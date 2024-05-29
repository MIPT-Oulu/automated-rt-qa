# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:32:35 2024

@author: rytkysan
"""
from pylinac import DRGS, DRMLC, WinstonLutz, ACRMRILarge
import numpy as np
import os
import logging
from pathlib import Path
from time import time

from pylinac import image

from qa_analysis.utilities import wait_user_close, move_file, save_excel


def drgs_test(mlc, open_im, tol=1.5, savepath=None, pdf=False, plot=False, precision=5,
              segment_size=None, roi=None, rep_dir='T2-T3 reports'):
    """
    Dose-rate and Gantry speed tests (T2 tests).
    Pylinac should automatically identify open beam and MLC images.
    Parameters
    ----------
    mlc : TYPE
        DESCRIPTION.
    open_im : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        Test tolerance (%). The default is 1.5.
    savepath : TYPE, optional
        DESCRIPTION. The default is None.
    pdf : bool, optional
        Toggle saving pdf report. The default is False.
    plot : bool, optional
        Toggle plotting the results image. The default is False.
    precision : int, optional
        Precision of measurement result. The default is 5 digits.
    segment_size: tuple, optional
        Sets the size for analysis segments in mm.
    roi: dict, optional
        Sets the offset positions and names for analysis segments.
    Returns
    -------
    dict
        Results of the analysis, rounded to precision digits.

    """
    
    # Test logger
    logger_t = logging.getLogger('qa.test')
    logger_t.info(f"Running DRGS test for {Path(mlc.path).name}")
    
    drgs = DRGS(image_paths=[open_im.path, mlc.path])
    
    if segment_size is not None and roi is not None:
        drgs.analyze(tolerance=tol, segment_size_mm=segment_size, roi_config=roi)
    elif roi is not None:
        drgs.analyze(tolerance=tol, roi_config=roi)
    else:
        drgs.analyze(tolerance=tol)
    
    # Plot figures
    if plot:
        drgs.plot_analyzed_image()
    
    
    # Save results
    if pdf and savepath is not None:
        report_name =f'{mlc.metadata.PatientID}_{mlc.metadata.SeriesDate}_{mlc.metadata.SeriesTime}_t3.pdf'
        (savepath / rep_dir).mkdir(exist_ok=True)  # Make reports directory
        path = str(savepath / rep_dir / report_name)
        if wait_user_close(path):
            drgs.publish_pdf(path, notes=[f'Device: {mlc.metadata.StationName}', f'Operator: {mlc.metadata.OperatorsName}'])
        
    res = drgs.results_data(as_dict=True)
    
    # Round to given precision
    res['max_deviation_percent'] = np.round(res['max_deviation_percent'], precision)
    res['abs_mean_deviation'] = np.round(res['abs_mean_deviation'], precision)
        
    return res
   
     
def drmlc_test(mlc, open_im, tol=1.5, savepath=None, pdf=False, plot=False, precision=5,
               segment_size=None, roi=None, rep_dir='T2-T3 reports'):
    """
    Dose-rate and MLC speed tests (T3 tests).
    Pylinac should automatically identify open beam and MLC images.

    Parameters
    ----------
    mlc : TYPE
        DESCRIPTION.
    open_im : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        Test tolerance (%). The default is 1.5.
    savepath : TYPE, optional
        DESCRIPTION. The default is None.
    pdf : bool, optional
        Toggle saving pdf report. The default is False.
    plot : bool, optional
        Toggle plotting the results image. The default is False.
    precision : int, optional
        Precision of measurement result. The default is 5 digits.
    segment_size: tuple, optional
        Sets the size for analysis segments in mm.
    roi: dict, optional
        Sets the offset positions and names for analysis segments.

    Returns
    -------
    dict
        Results of the analysis, rounded to precision digits.

    """
    # Test logger
    logger_t = logging.getLogger('qa.test')
    logger_t.info(f"Running DRMLC test for {Path(mlc.path).name}")
    
    drmlc = DRMLC(image_paths=[open_im.path, mlc.path])
    
    if segment_size is not None and roi is not None:
        drmlc.analyze(tolerance=tol, segment_size_mm=segment_size, roi_config=roi)
    elif roi is not None:
        drmlc.analyze(tolerance=tol, roi_config=roi)
    else:
        drmlc.analyze(tolerance=tol)
    
    # Plot figures
    if plot:
        drmlc.plot_analyzed_image()
    
    # Save results
    if pdf and savepath is not None:
        report_name = f'{mlc.metadata.PatientID}_{mlc.metadata.SeriesDate}_{mlc.metadata.SeriesTime}_t3.pdf'
        (savepath / rep_dir).mkdir(exist_ok=True)  # Make reports directory
        path = str(savepath / rep_dir / report_name)
        if wait_user_close(path):
            drmlc.publish_pdf(path, notes=[f'Device: {mlc.metadata.StationName}', f'Operator: {mlc.metadata.OperatorsName}'])
        
    res = drmlc.results_data(as_dict=True)
        
    # Round to given precision
    res['max_deviation_percent'] = np.round(res['max_deviation_percent'], precision)
    res['abs_mean_deviation'] = np.round(res['abs_mean_deviation'], precision)

    return res


def catphan_analysis(im, args, pdf=True, plot=False, tolerances=dict(), 
                     rep_dir='Catphan reports', timeout=5):
    """
    

    Parameters
    ----------
    im : TYPE
        DESCRIPTION.
    args : TYPE
        DESCRIPTION.
    pdf : TYPE, optional
        DESCRIPTION. The default is True.
    plot : TYPE, optional
        DESCRIPTION. The default is False.
    tolerances : TYPE, optional
        DESCRIPTION. The default is dict().
    rep_dir : TYPE, optional
        DESCRIPTION. The default is 'Catphan reports'.
    timeout : TYPE, optional
        Timeout for running the Catphan analysis. The default is 5 minutes.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    # Test logger
    logger_t = logging.getLogger('qa.test')
    logger_t.info(f"Running Catphan analysis for {Path(im.path).name}")
    
    
    analysis_path = os.path.dirname(im.path)
    start = time()
    
    while len(os.listdir(analysis_path)) > 0 and start - time() < timeout * 60:
        # Run the analysis for Catphan model assigned in args
        cbct = args.catphan_model(analysis_path)
       
        # Use the test tolerances from constants.py
        cbct.analyze(**tolerances)
        
        res = cbct.results_data(as_dict=True)
        res['mtf'] = cbct.ctp528.mtf.mtfs
        
        # Update DICOM metadata
        im = cbct.dicom_stack[0]
        
        # Add tolerances to results
        if len(tolerances) > 0:
            res['ctp404']['scaling_tolerance'] = tolerances['scaling_tolerance']
            res['ctp404']['thickness_tolerance'] = tolerances['thickness_tolerance']
            res['ctp404']['low_contrast_tolerance'] = tolerances['low_contrast_tolerance']
        else:
            # Pylinac defaults
            res['ctp404']['scaling_tolerance'] = 1
            res['ctp404']['thickness_tolerance'] = 0.2
            res['ctp404']['low_contrast_tolerance'] = 1
        
        # Plot figures
        if plot:
            cbct.plot_analyzed_image()
            
        # Save results
        if pdf:
            report_name = f'{im.metadata.PatientID}_{im.metadata.SeriesDate}_{im.metadata.SeriesTime}_Catphan.pdf'
            (args.save_path / rep_dir).mkdir(exist_ok=True)  # Make reports directory
            path = str(args.save_path / rep_dir / report_name)
            if wait_user_close(path):
                cbct.publish_pdf(path, notes=[f'Device: {im.metadata.StationName}', f'Operator: {im.metadata.OperatorsName}'])
        
        # Save Catphan analysis to Excel file
        save_excel(im, res, save_path=args.save_path, test='Catphan')
                
        # Move analyzed files to the processed folder, create subfolder by modality
        modality = 'Catphan'
        # Assume that there is one folder for patient name/ID
        parent_folder = Path(im.path).parent.parent.stem
        dicom_stack = image.DicomImageStack(analysis_path)
        for img in dicom_stack:        
            img.metadata[0x0018, 0x5100].value = 'HFS'
            # Replace the data folder in image path with processed
            if parent_folder == modality:
                processed_path = img.path.replace(args.data_path.stem, f'{args.processed_path.stem}' )
            else:
                processed_path = img.path.replace(args.data_path.stem, f'{args.processed_path.stem}/{modality}' )
            # Move the file
            move_file(img.path, processed_path)
            
    return res


def acr_analysis(im, args, pdf=True, plot=False, rep_dir='ACR reports', timeout=5):
    # Test logger
    logger_t = logging.getLogger('qa.test')
    logger_t.info(f"Running ACR analysis for {Path(im.path).name}")
    
    analysis_path = os.path.dirname(im.path)
    start = time()
    
    while len(os.listdir(analysis_path)) > 0 and start - time() < timeout * 60:
    
        dicom_stack = image.DicomImageStack(analysis_path, min_number=3)
        for img in dicom_stack:
            # Update field strength to Dicom metadata
            img.metadata.MagneticFieldStrength = args.field_strength
            # Save the image with new metadata
            img.save(img.path)    
    
        # Run the analysis for MR images of the ACR phantom
        acr = ACRMRILarge(os.path.dirname(im.path))    
        acr.analyze()
        
        # Plot figures
        if plot:
            acr.plot_analyzed_image()
            
        # Save results
        if pdf:
            report_name = f'{im.metadata.PatientID}_{im.metadata.SeriesDate}_{im.metadata.SeriesTime}_ACR.pdf'
            (args.save_path / rep_dir).mkdir(exist_ok=True)  # Make reports directory
            path = str(args.save_path / rep_dir / report_name)
            if wait_user_close(path):
                acr.publish_pdf(path)
                
        # Move analyzed files to the processed folder, create subfolder by modality
        modality = 'ACR'
        # Assume that there is one folder for patient name/ID
        parent_folder = Path(im.path).parent.parent.stem
        dicom_stack = image.DicomImageStack(analysis_path, min_number=3)
        for img in dicom_stack:        
            
            # Replace the data folder in image path with processed
            if parent_folder == modality:
                processed_path = img.path.replace(args.data_path.stem, f'{args.processed_path.stem}' )
            else:
                processed_path = img.path.replace(args.data_path.stem, f'{args.processed_path.stem}/{modality}' )
            
            # Move the file
            move_file(img.path, processed_path)
        
    return acr.results_data(as_dict=True)


def winston_analysis(im, args, pdf=True, plot=False, rep_dir='Winston-Lutz reports'):
        
    # Run the analysis for given image parent folder
    wl = WinstonLutz(os.path.dirname(im.path))
    
    # Test logger
    logger_t = logging.getLogger('qa.test')
    logger_t.info(f"Running Winston-Lutz analysis for {Path(im.path).name}")
    
    # Use the test tolerances from constants.py
    wl.analyze(bb_size_mm=args.bb_size_mm)
    
    # Plot figures
    if plot:
        wl.plot_summary()
        
    # Save results
    if pdf:
        report_name = f'{im.metadata.PatientID}_{im.metadata.StationName}_{im.metadata.SeriesDate}_{im.metadata.SeriesTime}_Winston_Lutz.pdf'
        (args.save_path / rep_dir).mkdir(exist_ok=True)  # Make reports directory
        path = str(args.save_path / rep_dir / report_name)
        if wait_user_close(path):
            wl.publish_pdf(path, notes=[f'Device: {im.metadata.StationName}', f'Operator: {im.metadata.OperatorsName}'])
        
    # Move analyzed files to the processed folder, create subfolder by modality
    modality = 'Winston-Lutz'
    # Assume that there is one folder for patient name/ID
    parent_folder = Path(im.path).parent.parent.stem
    # List files in the parent folder
    images = os.listdir(os.path.dirname(im.path))
    images.sort() 
    for img in images:     
        img = os.path.join(os.path.dirname(im.path), img)    
        
        # Replace the data folder in image path with processed
        if parent_folder == modality:
            processed_path = img.replace(args.data_path.stem, f'{args.processed_path.stem}' )
        else:
            processed_path = img.replace(args.data_path.stem, f'{args.processed_path.stem}/{modality}' )
        
        # Move the file
        move_file(img, processed_path)
        
    return wl.results_data(as_dict=True)