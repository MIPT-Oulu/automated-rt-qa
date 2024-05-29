# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:33:10 2024

@author: rytkysan
"""

import argparse
from pathlib import Path
from pylinac import (
    CatPhan503,
    CatPhan504,
    CatPhan600,
    CatPhan604
    )

from qa_analysis.constants import CustomCP504   
from qa_analysis.analysis import analyze_image
from qa_analysis.utilities import start_log

def main():
    # Input arguments and constants
    parser = argparse.ArgumentParser(
        description='Automated radiation therapy QA tests')
    parser.add_argument('--data_path', type=Path, default='Z:/Python/automated-rt-qa/data')
    parser.add_argument('--network_path', type=Path, default='share.txt', 
                        help='Path for a file with network drive details.')
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

    arg = parser.parse_args()
    
    # Set up logging for file and console
    start_log(arg.log_path)

    # Analysis script
    analyze_image(arg)
    
    
if __name__ == "__main__":   
    main()