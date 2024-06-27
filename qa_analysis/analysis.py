# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:32:35 2024

@author: rytkysan

Known issues: 
    - Identification of faulty image (cropped field)
    - Other data types than dicom
    - Pylinac 3.22 does not sort multiple series correctly. Issue raised:
        https://github.com/jrkerns/pylinac/issues/494
"""
from pylinac import image
from glob import glob
import numpy as np
from pathlib import Path
from os.path import dirname
import logging

from qa_analysis.tests import drgs_test, drmlc_test, catphan_analysis, winston_analysis, acr_analysis, normi_13_analysis
from qa_analysis.utilities import save_excel, move_file, remove_empty_dir, map_network_drive
from qa_analysis.constants import (
    T2_DR_ROI_HAL, T2_GS_ROI_HAL, T3_MLC_ROI_HAL, 
    DRGS_TOL, DRMLC_TOL, CATPHAN_CBCT_TOLERANCES, CATPHAN_TOLERANCES
    )
    

def analyze_image(arg):
    """
    Main analysis pipeline.
    
    Identifies T2-T3 test based on RT image label and Meterset exposure.
    With Halcyon DR test, exposure is ~60MU.
    
    Parameters
    ----------
    arg : TYPE
        Input arguments.

    Returns
    -------
    None.

    """
    # Analysis logger
    logger_a = logging.getLogger('qa.analysis')
    
    # Map network drive with correct login details
    if arg.network_path is not None:
        map_network_drive(arg.network_path)
    
    # List dicom files in data path
    images = []
    for file_type in arg.file_types:
        images += glob(str(arg.data_path / f'**/*{file_type}'), recursive=True)
    images.sort()
    
    # Check for empty directory
    if len(images) == 0:
        logger_a.info('No files in the analysis folder!')
        return
    
    # Load dicom images
    dcm_images = []
    dates = []  # Series date
    times = []  # Series time
    patients = []  # Patient ID (device name)
    files = []  # File names
    for i in range(len(images)):
        try:
            im = image.LinacDicomImage(images[i])
            dcm_images.append(im)
        except AttributeError:
            logger_a.info(f'Unable to open image {images[i]}')
            continue
        dates.append(im.metadata[0x0008, 0x0021].value)
        times.append(im.metadata[0x0008, 0x0031].value)
        patients.append(im.metadata[0x0010, 0x0020].value)
        files.append(Path(images[i]).stem)
    
    # Sort images by Series Date
    dates, times, dcm_images, patients, files = zip(
        *sorted(zip(dates, times, dcm_images, patients, files), 
                key=lambda x: x[0]))
    
    # Pick unique measurement dates, linac (and test?)
    unique_dates = np.unique(dates)
        
    # Loop for measurement dates
    for date in unique_dates:

        # Indexes of the measurement date
        date_id = np.where(np.array(dates) == date)[0]
        
        # Pick unique patients on the measurement date
        unique_patients = np.unique(np.array(patients)[date_id])
        
        # Loop for all patients
        for patient in unique_patients:            
            # Indexes of the same patient ID in the measurement date
            pat_id = np.where(np.array(patients) == patient)[0]
            # Combine the date and patient id
            pat_id = list(set(pat_id) & set(date_id))
                        
            try:
                # Loop through images of the patient in measurement date
                # Find the relevant images for each test
                test_images = {
                    't2_dr_segment_size': None,
                    't2_gs_segment_size': None,
                    't3_segment_size': None,
                    't2_dr_roi': None,
                    't2_gs_roi': None,
                    't3_roi': None}
                for im_id in pat_id:
                    # Find if the image is from T2 or T3 test
                    test_images = detect_t2_t3_tests(dcm_images[im_id], test_images)

                    # Find Catphan images
                    test_images = detect_catphan_tests(dcm_images[im_id], test_images)

                    # Find ACR images
                    test_images = detect_acr_tests(dcm_images[im_id], test_images, arg)

                    # Find Winston-Lutz images
                    test_images = detect_winston_tests(dcm_images[im_id], test_images)

                    # Find CR images
                    test_images = detect_normi_13_tests(dcm_images[im_id], test_images)

                # T2/T3 test detected
                if 't3_mlc' in test_images:
                    # Run T2/T3 analysis
                    run_t2_t3_tests(test_images, arg)
                elif 'catphan' in test_images:
                    # Run Caphan analysis
                    results = catphan_analysis(test_images['catphan'], arg,
                                               tolerances=CATPHAN_TOLERANCES)
                elif 'catphan_linac' in test_images:
                    # Run Caphan analysis
                    results = catphan_analysis(test_images['catphan_linac'], arg,
                                               tolerances=CATPHAN_CBCT_TOLERANCES)
                elif 'acr' in test_images:
                    # Run ACR analysis
                    results = acr_analysis(test_images['acr'], arg)
                elif 'winston' in test_images:
                    # Run Winston-Lutz analysis
                    results = winston_analysis(test_images['winston'], arg)
                elif 'normi_13' in test_images:
                    results = normi_13_analysis(test_images['normi_13'], arg)
                # Here more tests could be run
                else:
                    logger_a.info(f'Test not implemented for patient {patient}, date {date}')

            # Missing dictionary data raises KeyError
            # ValueError when running Winston analysis with incorrect images
            except (KeyError, TypeError, ValueError,
                    AttributeError, ZeroDivisionError, IndexError) as e:
                logger_a.debug(f'Cannot analyse from measurement date {date} due to error {e}')

    # List dicom files remaining in data path
    images = glob(str(arg.data_path / '**/*.*'), recursive=True)
    images.sort()
    # Move files to the processed folder
    for im_id in range(len(images)):
        # Replace the data folder in image path with processed
        processed_path = images[im_id].replace(arg.data_path.stem, f'{arg.processed_path.stem}/Not_analyzed')
        # Move the file
        move_file(images[im_id], processed_path)
        
    # Check for empty directories in data path
    remove_empty_dir(arg.data_path)


def detect_t2_t3_tests(dcm_image, res_images):
    """
    Finds open-beam and MLC images for T2 and T3. 
    Sets custom ROI properties for Halcyon.

    Assumptions:
    T2: RT image label = 'MV_243'
    T3: RT image label = 'MV_32' (MV_190 open and MV_40 MLC with Halcyon)
    Open beam images have the "Curve label" -tag (= 'Field Edge (Open' or 'CIAO (OpenBeam)' with Halcyon)
    Halcyon T2DR: Exposure < 100 MU
    Halcyon T2GS: Jaw position = -140 mm


    Parameters
    ----------
    dcm_image : TYPE
        DESCRIPTION.
    res_images : TYPE
        DESCRIPTION.

    Returns
    -------
    res_images : TYPE
        DESCRIPTION.

    """
    # Check if RT image label exists in Dicom metadata
    if (0x3002, 0x0002) in dcm_image.metadata:
        rt_label = dcm_image.metadata[0x3002, 0x0002].value
    else:
        return res_images

    # T2 open beam
    if 'MV_243' in rt_label and (0x5000, 0x2500) in dcm_image.metadata:
        # Halcyon DR image has ~60MU exposure
        exposure = float(dcm_image.metadata[0x3002, 0x0030].value[0].MetersetExposure)
        if exposure < 100:
            res_images['t2_dr_open'] = dcm_image
        else:
            res_images['t2_open'] = dcm_image
    # T2 MLC
    elif 'MV_243' in rt_label and not (0x5000, 0x2500) in dcm_image.metadata:
        # Halcyon DR image has ~60MU exposure
        exposure = float(dcm_image.metadata[0x3002, 0x0030].value[0].MetersetExposure)
        if exposure < 100:
            res_images['t2_dr_mlc'] = dcm_image
            # Custom ROI for Halcyon T2DR  
            res_images['t2_dr_roi'] = T2_DR_ROI_HAL
        else:
            res_images['t2_mlc'] = dcm_image
            
            # Halcyon T2 tests have 14cm x 14cm collimation
            jaw_pos = int(dcm_image.metadata[0x3002, 0x0030].value[0][0x300a, 0x00b6].value[0].LeafJawPositions[0])
            # Custom ROI for Halcyon T2GS
            if jaw_pos == -140: 
                res_images['t2_gs_roi'] = T2_GS_ROI_HAL
    # T3 Open beam
    elif ('MV_32' in rt_label or 'MV_190' in rt_label) and (0x5000, 0x2500) in dcm_image.metadata:
        res_images['t3_open'] = dcm_image
        
        # Custom ROI for Halcyon T3MLC
        if 'MV_190' in rt_label:                                                     
            res_images['t3_roi'] = T3_MLC_ROI_HAL
    # T3 MLC
    elif ('MV_32' in rt_label or 'MV_40' in rt_label) and not (0x5000, 0x2500) in dcm_image.metadata:
        res_images['t3_mlc'] = dcm_image

    return res_images


def detect_catphan_tests(dcm_image, res_images):

    # Modality should be CT
    if dcm_image.metadata[0x0008, 0x0060].value != 'CT':
        return res_images
    
    # The orientation should be head first supine
    if dcm_image.metadata[0x0018, 0x5100].value == 'FFS':
        # Invert Z-axis, and change to HFS
        dcm_image.metadata[0x0018, 0x5100].value = 'HFS'
        dcm_image.metadata[0x0020, 0x0032].value[2] = str(-float(dcm_image.metadata[0x0020, 0x0032].value[2]))
        # Save inverted image
        dcm_image.save(dcm_image.path)
    
    # Linac CBCT
    if (0x0008, 0x114a) in dcm_image.metadata:
        
        # For Linac CBCT, SOP UID should be RT Plan storage
        ref_inst = 'RT Plan or RT Ion Plan or Radiation Set to be verified'
        catphan_test1 = dcm_image.metadata[0x0008, 0x114a].value[0][0x0008, 0x1150].repval == 'RT Plan Storage'
        catphan_test2 = dcm_image.metadata[0x0008, 0x114a].value[0][0x0040, 0xa170].value[0][0x0008, 0x0104].value == ref_inst
        
        if catphan_test1 or catphan_test2:
            
            
            # Save the analysis image if does not exist already
            if not 'catphan_linac' in res_images:
                res_images['catphan_linac'] = dcm_image
    
    # Diagnostic CT
    elif (0x0008, 0x1140) in dcm_image.metadata and not 'catphan' in res_images:
        # For diagnostic CT, SOP UID [0x0008, 0x1140][0x0008, 0x1150] should be CT Image Storage
        if dcm_image.metadata[0x0008, 0x1140].value[0][0x0008, 0x1150].repval == 'CT Image Storage':
            res_images['catphan'] = dcm_image
        
    return res_images


def detect_acr_tests(dcm_image, res_images, args):
    """
    Saves the first MR image found to given dict.
    Assumes that MR images are from ACR phantom 
    and that the images from same series are in one folder.

    Parameters
    ----------
    dcm_image : LinacDicomImage
        DESCRIPTION.
    res_images : dict
        DESCRIPTION.
    args : TYPE
        Input arguments.

    Returns
    -------
    res_images : dict
        DESCRIPTION.

    """

    # Modality should be MRI
    if dcm_image.metadata[0x0008, 0x0060].value != 'MR':
        return res_images

    if 'acr' not in res_images:
        # Update field strength to Dicom metadata
        dcm_image.metadata.MagneticFieldStrength = args.field_strength
        # Save the image with new metadata
        dcm_image.save(dcm_image.path)
        
        res_images['acr'] = dcm_image
    # ACR image from same series
    elif dirname(dcm_image.path) == dirname(res_images['acr'].path):
        pass

    return res_images


def detect_winston_tests(dcm_image, res_images):

    # Patient name could be added as a filter    

    # Modality should be RTIMAGE
    if dcm_image.metadata[0x0008, 0x0060].value != 'RTIMAGE':
        return res_images

    # First WL test image
    if not 'winston' in res_images:
        res_images['winston'] = dcm_image
    # Winston-Lutz image from same series
    elif dirname(dcm_image.path) == dirname(res_images['winston'].path):
        pass
        
    return res_images


def detect_normi_13_tests(dcm_image, res_images):
    # Modality should be CR
    modality = dcm_image.metadata[0x0008, 0x0060].value
    if (modality == 'DX' or modality == 'CR') and 'cr' not in res_images:
        res_images['normi_13'] = dcm_image

    return res_images


def run_t2_t3_tests(test, args):   
    res = []
    
    # Dose-rate & gantry speed test (T2)
    t2 = drgs_test(test['t2_mlc'], test['t2_open'], tol=DRGS_TOL, 
              savepath=args.save_path, pdf=args.pdf, plot=args.plot,
              segment_size=test['t2_gs_segment_size'], roi=test['t2_gs_roi'])
    res.append(t2)
    
    # mlc speed test (T3)
    t3 = drmlc_test(test['t3_mlc'], test['t3_open'], tol=DRMLC_TOL, 
              savepath=args.save_path, pdf=args.pdf, plot=args.plot,
              segment_size=test['t3_segment_size'], roi=test['t3_roi'])
    res.append(t3)
    
    # Dose rate test for Halcyon
    if 't2_dr_open' and 't2_dr_mlc' in test:
        t2_dr = drgs_test(test['t2_dr_mlc'], test['t2_dr_open'], tol=DRGS_TOL, 
                  savepath=args.save_path, pdf=args.pdf, plot=args.plot,
                  segment_size=test['t2_dr_segment_size'], roi=test['t2_dr_roi'])
        res.append(t2_dr)
    
    # Save results as a row in Excel file
    save_excel(test['t2_mlc'], res, save_path=args.save_path, test='T2-T3')
    
    # Move analyzed files to the processed folder, create subfolder by modality
    modality = 'T2-T3'
    # Assume that there is one folder for patient name/ID
    parent_folder = Path(test['t2_mlc'].path).parent.parent.stem  
    # Possible test images
    t2t3_images = ['t2_mlc', 't2_open', 't3_mlc', 't3_open', 't2_dr_mlc', 't2_dr_open']
    for key, im in test.items():
        if key in t2t3_images:            
            # Replace the data folder in image path with processed
            if parent_folder == modality:
                processed_path = im.path.replace(args.data_path.stem, f'{args.processed_path.stem}' )
            else:
                processed_path = im.path.replace(args.data_path.stem, f'{args.processed_path.stem}/{modality}' )
            # Move the file
            move_file(im.path, processed_path)
