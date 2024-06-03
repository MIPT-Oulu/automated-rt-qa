# Automated quality control for radiation therapy 

(c) Santeri Rytky, Oulu University Hospital, 2024

## Background
Multiple quality assurance tests are done periodically at radiation therapy clinics. 
The analysis of the measurements is often time-consuming and subject to errors from manual processing.
This repository acts as a framework that automatically analyses images from a given data folder.
The main dependencies are the Pylinac package for calculating the analyses 
and the Watchdog package for monitoring the data folder.

## Requirements
There are small modifications to the Pylinac 3.22 package that should be done.
These are likely fixed in future updates.

Fix for reading multiple CT series ([Github issue](https://github.com/jrkerns/pylinac/issues/494))

### Installation
- [Anaconda installation](https://docs.anaconda.com/anaconda/install/) 
```
git clone https://github.com/MIPT-Oulu/automated-rt-qa.git
cd automated-rt-qa
```
Either create a conda environment for running the code
```
conda create --name rt-env --file requirements.txt
```
or install the packages using pip
```
pip install -r requirements.txt
```

Populate the input arguments from `main.py` and `main_offline.py`. 
If you need to map a network drive for the data folder, include a `share.txt` file.
The file should include the drive letter, shared path, username and password. Example:
```
Z:
\\network\drive
user
pwd
```

## Usage
The automated analysis is started using `main.py`.
This could be automated for example with task scheduler in Windows systems.
Other option is to run the analysis using `main_offline.py`, which runs the analysis pipeline once.

## Features

### Automated QA pipeline
The code lists the files of listed data types (default: .dcm, .tif, .tiff) in the folder, and 
runs the analysis for each measurement date (Series date) and Patient ID (test patient).
The code detects the available test types and runs all the different measurements found.
The results are saved either as a pdf report and/or a row in an Excel file.

### Available tests
- VMAT (T2/T3)
- Catphan analysis
- ACR phantom analysis (for MRI)
- Winston-Lutz

### Logging
Different events during the analysis pipeline are logged in the repository root.

## Assumptions
The pipeline is tested for Varian Truebeam and Halcyon accelerators.

### VMAT (T2/T3) analysis
The pipeline assumes that T2 images have RT label `MV_243`. 
In T3 tests the label is `MV_32` or `MV_190`.
Open beam images have the 5000, 2500 (Curve label) DICOM tag.

For Halcyon, the T2 test is split into two separate tests. 
The Halcyon dose-rate test is detected from a lower exposure (less than 100 MU).
Halcyon tests are conducted with custom regions of interest (`constants.py`).

### Catphan analysis
The DICOM modality attribute is `CT`.
CBCT images have either `RT Plan Storage` or `RT Plan or RT Ion Plan or Radiation Set to be verified` identifier.
Diagnostic CT images are detected based on `CT Image Storage` identifier.
Phantom is oriented head first supine (HFS). If feet first supine is planned, the code inverts the stack Z-axis.

### ACR analysis
Currently only passed for MR images. The modality should be `MR`.

### Winston-Lutz analysis
Modality should be `RTIMAGE`. The test is tried for images that are not labeled as other test types.

## License
This software is distributed under the MIT License.
