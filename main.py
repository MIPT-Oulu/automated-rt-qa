# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:33:10 2024

@author: rytkysan
"""

import argparse
import logging
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
from pathlib import Path
from time import sleep
import os
from tqdm import tqdm
from pylinac import (
    CatPhan503,
    CatPhan504,
    CatPhan600,
    CatPhan604
    )

from qa_analysis.constants import CustomCP504    
from qa_analysis.analysis import analyze_image
from qa_analysis.utilities import map_network_drive, start_log


def main():
    # Input arguments and constants
    parser = argparse.ArgumentParser(
        description='Automated radiation therapy QA tests')
    parser.add_argument('--data_path', type=Path, default='Z:/Python/automated-rt-qa/data')
    parser.add_argument('--network_path', type=Path, default='share.txt')
    parser.add_argument('--processed_path', type=Path, default='Z:/Python/automated-rt-qa/processed')
    parser.add_argument('--save_path', type=Path, default='Z:/Python/automated-rt-qa/results')
    parser.add_argument('--log_path', type=Path, default='automated_qa.log', help='File for saving event logs.')
    parser.add_argument('--file_types', type=tuple, default=('.dcm', '.tiff', '.tif'), help='File types listed for analysis.')
    parser.add_argument('--catphan_model', default=CustomCP504, 
                        choices=[CatPhan503, CatPhan504, CatPhan600, CatPhan604, CustomCP504], 
                        help='Catphan phantom model')
    parser.add_argument('--field_strength', type=float, default=3.0, help='MRI field strength for ACR phantom images.')
    parser.add_argument('--bb_size_mm', type=float, default=6.0, 
                        help='Size of the ball-bearing phantom for Winston-Lutz test.')
    parser.add_argument('--pdf', type=bool, default=False, help='Option for saving a pdf results file.')
    parser.add_argument('--plot', type=bool, default=False, help='Option for plotting results images.')    
    parser.add_argument('--wait_time', type=int, default=30, help='Waiting time (s) after a file is found. Allows user to finish file transfers.')
    parser.add_argument('--monitor_time', type=int, default=5, help='Waiting time (s) for checking if file structure has changed.')
    
    # Use a global variable for arguments to allow updating them outside the function
    global arg
    arg = parser.parse_args()
    
    # Map network drive with correct password
    if arg.network_path is not None:
        map_network_drive(arg.network_path)
    
    # --- Watch dog ---
    # Define Watchdog event handler
    automated_qa = AutomatedQA(patterns="*",
                               ignore_patterns="",
                               ignore_directories=False,
                               case_sensitive=True)
    
    # Set up logging for file and console
    start_log(arg.log_path)
    
    # Watchdog observer to monitor data folder
    observer = Observer()
    observer.schedule(automated_qa, arg.data_path, recursive=True)
    observer.start()
    
    try:
        while True:
            sleep(arg.monitor_time)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

class AutomatedQA(PatternMatchingEventHandler):
    def on_created(self, event):
        # Skip the queue when directory is empty
        if len(os.listdir(str(arg.data_path))) == 0:
            return
            
        # Log the first file found
        logging.info(f"{event.src_path} found. Running analysis...")
        
        # Wait to allow finishing file transfers
        division = 100
        for i in tqdm(range(division), desc=f'Wait {arg.wait_time}s to allow finishing file transfers'):
            sleep(arg.wait_time / division)
        
        # Run the analysis
        analyze_image(arg)
        
    def on_deleted(self, event):
        # Log only processing of directories
        if event.is_directory:
            logging.info(f"{event.src_path} is deleted.")
    
    def on_modified(self, event):
        # Log only processing of directories
        if event.is_directory:
            logging.info(f"{event.src_path} is modified.")
    
    def on_moved(self, event):
        # Log only processing of directories
        if event.is_directory:
            logging.info(f"Moved from {event.src_path} to {event.dest_path}.")
        

if __name__ == "__main__":   
    main()