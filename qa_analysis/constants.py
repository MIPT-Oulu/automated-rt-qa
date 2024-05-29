# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:30:18 2024

@author: rytkysan
"""

from pylinac.ct import CatPhan504, CTP404CP504, CTP486, CTP528CP504, CTP515

# ROI sizes for Halcyon T2/T3 tests

T2_DR_ROI_HAL = {
    '-6 cm': {'offset_mm': -60},
    '-4 cm': {'offset_mm': -40},
    '-2 cm': {'offset_mm': -20},
    '0 cm': {'offset_mm': 0},
    '+2 cm': {'offset_mm': 20},
    '+4 cm': {'offset_mm': 40},
    '+6 cm': {'offset_mm': 60}}

T2_GS_ROI_HAL = {
    '-12 cm': {'offset_mm': -120},
    '-8 cm': {'offset_mm': -80},
    '-4 cm': {'offset_mm': -40},
    '0 cm': {'offset_mm': 0},
    '+4 cm': {'offset_mm': 40},
    '+8 cm': {'offset_mm': 80},
    '+12 cm': {'offset_mm': 120},}

T3_MLC_ROI_HAL = {
    '-11.2 cm': {'offset_mm': -112},
    '-5.6 cm': {'offset_mm': -56},                                
    '0 cm': {'offset_mm': 0},
    '+5.6 cm': {'offset_mm': 56},
    '+11.2 cm': {'offset_mm': 112}}

# Tolerances for Catphan analysis

CATPHAN_TOLERANCES = {
    'hu_tolerance': 10,
    'scaling_tolerance': 0.5,
    'thickness_tolerance': 0.1,
    'low_contrast_tolerance': 5,
    'expected_hu_values': {'Air': -983, 
                           'PMP': -181, 
                           'LDPE': -92, 
                           'Poly': -37, 
                           'Acrylic': 122, 
                           'Delrin': 343, 
                           'Teflon': 936}}

CATPHAN_CBCT_TOLERANCES = {
    'hu_tolerance': 20,
    'scaling_tolerance': 0.5,
    'thickness_tolerance': 0.1,
    'low_contrast_tolerance': 5,
    'expected_hu_values': {'Air': -983, 
                           'PMP': -181, 
                           'LDPE': -92, 
                           'Poly': -37, 
                           'Acrylic': 122, 
                           'Delrin': 343, 
                           'Teflon': 936}}

# Tolerance for DRGS (T2) test (% of max deviation)
DRGS_TOL = 1.5

# Tolerance for DRMLC (T3) test (% of max deviation)
DRMLC_TOL = 1.5

# Custom linearity module (smaller diameter for ROIs)

AIR = -1000
PMP = -196
LDPE = -104
POLY = -47
ACRYLIC = 115
DELRIN = 365
TEFLON = 1000
BONE_20 = 237
BONE_50 = 725
WATER = 0

class CustomCTP404(CTP404CP504):
    roi_dist_mm = 58.7  # Default value
    roi_radius_mm = 4  # Smaller diameter
    roi_settings = {
        "Air": {
            "value": AIR,
            "angle": -90,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "PMP": {
            "value": PMP,
            "angle": -120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "LDPE": {
            "value": LDPE,
            "angle": 180,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Poly": {
            "value": POLY,
            "angle": 120,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Acrylic": {
            "value": ACRYLIC,
            "angle": 60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Delrin": {
            "value": DELRIN,
            "angle": 0,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
        "Teflon": {
            "value": TEFLON,
            "angle": -60,
            "distance": roi_dist_mm,
            "radius": roi_radius_mm,
        },
    }
    background_roi_settings = {
        "1": {"angle": -30, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "2": {"angle": -150, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "3": {"angle": -210, "distance": roi_dist_mm, "radius": roi_radius_mm},
        "4": {"angle": 30, "distance": roi_dist_mm, "radius": roi_radius_mm},
    }


# then, pass to the CatPhan model
class CustomCP504(CatPhan504):
    modules = {
        CustomCTP404: {"offset": 0},
        CTP486: {"offset": -65},
        CTP528CP504: {"offset": 30},
        CTP515: {"offset": -30},
    }