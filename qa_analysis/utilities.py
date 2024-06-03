# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
from pathlib import Path
from time import sleep, time
from subprocess import run


def start_log(path):
    # Log folder
    path.parent.mkdir(exist_ok=True)
    
    # Logging parameters
    month = datetime.today().strftime('_%Y_%m')  # Name log files using current year and month
    logging.basicConfig(filename=f'{path.parent}/{path.stem}{month}{path.suffix}', 
                        level=logging.DEBUG, 
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filemode='a'  # a = append to existing file, w = overwrite
                        )
    
    # Define a console Handler which writes INFO messages or higher
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # Set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # Tell the handler to use this format
    console.setFormatter(formatter)
    # Add the handler to the root logger
    logging.getLogger('').addHandler(console)


def map_network_drive(path):
    """
    Map network drive with address and credentials given in a file.

    Parameters
    ----------
    path : TYPE
        File path.

    Returns
    -------
    None.

    """
    # Map network drive with correct password
    with open(str(path), 'r') as f:
        share = f.readlines()
    
    # Return if path exists
    if os.path.exists(f'{share[0]}/'):
        return
    
    # Otherwise, map drive
    share_command=fr'net use {share[0]} {share[1]} /user:{share[2]} {share[3]}'.replace('\n', ' ')
    run(share_command, shell=True)

def save_excel(dicom_im, res, save_path, test='T2-T3', prec=5):
    """
    

    Parameters
    ----------
    dicom_im : TYPE
        Dicom image analyzed (for extracting test metadata).
    res : dict
        T2, T3 and T2_dr results
    save_path : TYPE
        DESCRIPTION.
    test : TYPE, optional
        QA test to be saved. The default is 'T2-T3'. 'Catphan' is also available.
    prec : int, optional
        Numeric precision for floating point results. The default is 5.

    Raises
    ------
    NotImplementedError
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # Date, time and Patient ID
    date = dicom_im.metadata[0x0008, 0x0021].value
    time = dicom_im.metadata[0x0008, 0x0031].value
    patient = dicom_im.metadata[0x0010, 0x0020].value
    
    # Results column headers
    if test == 'T2-T3':
        # Round the results to five digits
        t2_res, t3_res = res[0], res[1]     
        
        # Columns should be: Series date, Series time, Patient ID, Test, Pass/Fail, Max_deviation, Average_deviation
        cols = ['Series date', 
                'Series time', 
                'T2_Pass/Fail', 'T2_Max_deviation',
                'T3_Pass/Fail', 'T3_Max_deviation',
                'T2 DR GS (Avg)', 'T3 MLC SPEED (Avg)',]

        # Row of test results, in Excel-friendly format
        results_data = [f'{date[6:8]}.{date[4:6]}.{date[:4]}',
                        f'{time[:2]}:{time[2:4]}:{time[4:6]}',
                        t2_res['passed'], t2_res['max_deviation_percent'],
                        t3_res['passed'], t3_res['max_deviation_percent'], 
                        t2_res['abs_mean_deviation'], t3_res['abs_mean_deviation']]  
        
        # Halcyon data
        if len(res) == 3:
            t2dr_res = res[2]
            
            # Column names for Halcyon
            cols = ['Series date', 'Series time', 
                    'T2DR_Pass/Fail', 'T2DR_Max_deviation', 
                    'T2GS_Pass/Fail', 'T2GS_Max_deviation',
                    'T3_Pass/Fail', 'T3_Max_deviation',                    
                    'T2 DR (Avg)', 'T2 GS (Avg)', 'T3 LS (Avg)']
            
            # Row of test results, in Excel-friendly format
            results_data = [f'{date[6:8]}.{date[4:6]}.{date[:4]}', f'{time[:2]}:{time[2:4]}:{time[4:6]}',
                            t2dr_res['passed'], t2dr_res['max_deviation_percent'], 
                            t2_res['passed'], t2_res['max_deviation_percent'], 
                            t3_res['passed'], t3_res['max_deviation_percent'],                            
                            t2dr_res['abs_mean_deviation'], t2_res['abs_mean_deviation'], t3_res['abs_mean_deviation'],]  
            
    elif test == 'Catphan':
        tols = [res['ctp404']['hu_tolerance'], 
                res['ctp404']['scaling_tolerance'],
                res['ctp404']['thickness_tolerance'],
                res['ctp404']['low_contrast_tolerance']]
        
        # TODO what if constants are changed (tolerances)?
        lps, mtfs = zip(*res['mtf'].items()) 
        cols = ['Series date', 
                'Series time', 
                'Series description',
                'KVP',
                'mAs',
                'Filter type',
                'Convolution kernel',
                'CTDIvol',
                f'Linearity (HU, +/-{tols[0]}), Air',
                'PMP',
                'LDPE',
                'Polystyrene',
                'Acrylic',
                'Delrin',
                'Teflon',
                f'Average line distance (mm, < {tols[1]}mm error)',
                f'Slice thickness (mm, < {tols[2]}mm error)',                
                'Uniformity (HU, +/-{tols[0]}), Center',
                'Top',
                'Right',
                'Bottom',
                'Left',
                'Low contrast visibility',
                f'Low contrast ROIs seen (> {tols[3]})',
                'MTF 80%',
                'MTF 50%',
                'MTF 30%',
                f'MTF {lps[0]} lp/mm',
                f'MTF {lps[1]} lp/mm',
                f'MTF {lps[2]} lp/mm',
                f'MTF {lps[3]} lp/mm',
                f'MTF {lps[4]} lp/mm',
                f'MTF {lps[5]} lp/mm',
                f'MTF {lps[6]} lp/mm',
            ]
        
        series = ''
        ctdi = round(dicom_im.metadata[0x0018, 0x9345].value, prec) if (0x0018, 0x9345) in dicom_im.metadata else ''
        
        # Row of test results, in Excel-friendly format
        results_data = [f'{date[6:8]}.{date[4:6]}.{date[:4]}',
                        f'{time[:2]}:{time[2:4]}:{time[4:6]}',
                        series,
                        int(dicom_im.metadata[0x0018, 0x0060].value),
                        int(dicom_im.metadata[0x0018, 0x1152].value),
                        dicom_im.metadata[0x0018, 0x1160].value,
                        dicom_im.metadata[0x0018, 0x1210].value,
                        ctdi,
                        # Linearity, difference from reference
                        int(res['ctp404']['hu_rois']['Air']['difference']),
                        int(res['ctp404']['hu_rois']['PMP']['difference']),
                        int(res['ctp404']['hu_rois']['LDPE']['difference']),
                        int(res['ctp404']['hu_rois']['Poly']['difference']),
                        int(res['ctp404']['hu_rois']['Acrylic']['difference']),
                        int(res['ctp404']['hu_rois']['Delrin']['difference']),
                        int(res['ctp404']['hu_rois']['Teflon']['difference']),
                        # Geometry
                        round(res['ctp404']['avg_line_distance_mm'], prec), # Avg line distance
                        round(res['ctp404']['measured_slice_thickness_mm'], prec), # Slice thickness
                        # Uniformity
                        int(res['ctp486']['rois']['Center']['value']),
                        int(res['ctp486']['rois']['Top']['value']),                      
                        int(res['ctp486']['rois']['Right']['value']),
                        int(res['ctp486']['rois']['Bottom']['value']),                        
                        int(res['ctp486']['rois']['Left']['value']),
                        # Low contrast
                        round(res['ctp404']['low_contrast_visibility'], prec),  # Contrast of LDPE and Poly in the linearity module
                        res['ctp515']['num_rois_seen'],  # How many low-contrast ROIs detected
                        # MTF
                        round(res['ctp528']['mtf_lp_mm'][80], prec), # MTF 80%
                        round(res['ctp528']['mtf_lp_mm'][50], prec), # MTF 50% (Half-power frequency)
                        round(res['ctp528']['mtf_lp_mm'][30], prec), # MTF 30%
                        round(mtfs[0], prec), # MTF values for specific line pairs
                        round(mtfs[1], prec),
                        round(mtfs[2], prec),
                        round(mtfs[3], prec),
                        round(mtfs[4], prec),
                        round(mtfs[5], prec),
                        round(mtfs[6], prec),
                        ]
    else:
        raise NotImplementedError()
    
    # Compile results into dataframe
    results = pd.DataFrame(columns=cols)
    results.loc[0] = results_data
  
    
    # Add a new row to the excel file
    path_excel = str(save_path / f'Results_{patient}.xlsx')
    
    
        
    # Check if a results file exists
    if os.path.isfile(path_excel):
        
        # Check if the excel is opened
        if wait_user_close(path_excel):
            
            # Append to existing results file
            with pd.ExcelWriter(path_excel, mode='a', engine='openpyxl', if_sheet_exists='overlay') as writer:    
                
                # Create a new sheet to Excel
                if not test in writer.book:
                    results.to_excel(writer, sheet_name=test, index=None)
                
                # Find if the results row exists
                if find_matching_row(writer.book[test], results):
                    # Utility logger
                    logger_u = logging.getLogger('qa.utilities')
                    logger_u.info(f'Measurement date {date}, patient {patient}, test {test} already analyzed.')
                    return
                
                # Find the last row in the results table
                start = writer.book[test].max_row
                
                # Test for more parameters
                if len(results.columns) > writer.book[test].max_column:
                    for i, header in enumerate(results.columns):
                        writer.book[test].cell(1, i + 1).value = header
                    
                # Continue from the next row
                results.to_excel(writer, sheet_name=test, 
                                 header=None, index=None, startrow=start)
                
        # After timeout, skip saving the results excel
        else:        
            return
    
    # Create new results file
    else:            
        with pd.ExcelWriter(path_excel, engine='openpyxl') as writer:    
            results.to_excel(writer, sheet_name=test, index=None)


def move_file(src: str, dst: str, overwrite=True):
    """
    Moves the given file to a new destination. 
    Creates any folders that are needed for the destination.

    Parameters
    ----------
    src : str
        File to be moved.
    dst : str
        Destination for the file.
    Returns
    -------
    None.

    """
    # Utility logger
    logger_u = logging.getLogger('qa.utilities')
    
    # Create destination directory
    if not os.path.isdir(dst):
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
    
    # Move file
    # Check if file already exists
    if os.path.isfile(dst):
        # Possible to overwrite
        if overwrite:
            # Remove dst file and move src file to dst
            try:
                os.remove(dst)    
            except PermissionError:
                logger_u.debug(f'Unable to delete file {dst}')
                return
            os.rename(src, dst)
        else:
            logger_u.info(f'There already exists a file {dst}')
    else:
        # Move file
        os.rename(src, dst)
        
def remove_empty_directory(directory: Path):
    """
    Remove empty directories in a given path. Deprecated.

    Parameters
    ----------
    directory : Path
        Path for removing empty directories.

    Returns
    -------
    None.

    """
    for item in os.listdir(directory):
        # Full path
        item = str(directory / item)
        # Check for directory
        if os.path.isdir(item):
            # Remove empty directory
            if len(os.listdir(item)) == 0:
                os.rmdir(item)
                
def remove_empty_dir(pathlib_root_dir : Path):
    """
    Recursively remove empty directories in a given path.

    Parameters
    ----------
    directory : Path
        Path for removing empty directories.

    Returns
    -------
    None.

    """
    
    # List all directories recursively and sort them by path, longest first
    L = sorted(
        pathlib_root_dir.glob("**"),
        key=lambda p: len(str(p)),
        reverse=True,
    )
    # Do not remove the parent directory
    if pathlib_root_dir in L:
        L.remove(pathlib_root_dir)
    
    for pdir in L:
      try:
        pdir.rmdir()  # Remove directory if empty
      except OSError:
        continue  # Catch and continue if non-empty
        
    
def wait_user_close(path_file, retry_time=5, timeout=120):
    """
    Tries to open the file in a given path and tests if a user has 
    already opened the file. Returns True if able to open, False otherwise.

    Parameters
    ----------
    path_file : str
        Path to the file that the function tries to open. 
        If it is open, script retries after 5 seconds
    retry_time : int, optional
        Time to wait until retry for opening the file.
        The default is 5 seconds.
    timeout : int, optional
        Timeout for waiting the file to be closed. 
        The default is 120 minutes.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    # Utility logger
    logger_u = logging.getLogger('qa.utilities')
    
    if not os.path.isfile(path_file):
        return True
    
    s_time = time()
    
    # Wait until the Excel is closed
    while True:
        try:
            open(path_file, 'r+')
            break
        except PermissionError:                
            logger_u.info(path_file, ' is already opened. Waiting user to close...')
            sleep(retry_time)
            
            
            if time() - s_time > timeout * 60:
                logger_u.debug(f'Timeout of {timeout} minutes has passed. Returning...')
                return False
    return True

def find_matching_row(worksheet, compare_row):
    """
    Finds 

    Parameters
    ----------
    worksheet : TYPE
        DESCRIPTION.
    compare_row : TYPE
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    
    # Convert row to tuple
    compare_row = tuple(map(tuple, compare_row.values))[0]
    
    # Iterate through rows of the excel file to find a matching row
    for row in worksheet.iter_rows(values_only=True):
        if compare_row == row:
            return True
    return False

