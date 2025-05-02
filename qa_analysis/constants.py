# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:30:18 2024

@author: rytkysan
"""

from __future__ import annotations

from pylinac.ct import CatPhan504, CTP404CP504, CTP486, CTP528CP504, CTP515
from pylinac.cheese import CIRS062M, CheeseModule, CheeseResult
from pylinac.core.roi import HighContrastDiskROI
from pathlib import Path
import matplotlib.pyplot as plt
from typing import BinaryIO, Callable

from functools import cached_property
import numpy as np
from skimage.measure._regionprops import RegionProperties
from pylinac.ct import combine_surrounding_slices
from skimage import filters, draw, measure, segmentation
from scipy import ndimage
from pylinac.core.image import ImageLike
from pylinac.core import image
from pylinac.core.geometry import Point
from pylinac.core.profile import CollapsedCircleProfile
from pylinac.core.roi import RectangleROI
from pylinac.core.nps import plot_nps1d, noise_power_spectrum_2d, noise_power_spectrum_1d
from pylinac.vmat import VMATBase

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

# RapidArc dynamic test
RAD_ROI = {
    'A': {'offset_mm': -56.3, 'offset_y_mm': 34},
    'B': {'offset_mm': -56.3, 'offset_y_mm': -34},                                
    'C': {'offset_mm': 0, 'offset_y_mm': -66.5},
    'D': {'offset_mm': 56.3, 'offset_y_mm': 34},
    'E': {'offset_mm': 56.3, 'offset_y_mm': -34}}

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


class RAD(VMATBase):
    """Class representing RapidArc dynamic commissioning test. Will accept, analyze, and return the results."""

    _result_header = "RapidArc Dynamic synchronization test"
    _result_short_header = "RAD sync"
    default_roi_config = RAD_ROI
    
    def _calculate_segment_centers(self) -> list[Point]:
        """Construct the center points of the segments based on the field center and known x-offsets."""
        points = []
        _, open_prof = self._median_profiles(self.dmlc_image, self.open_image)
        x_field_center = round(open_prof.center_idx)
        for roi_data in self.roi_config.values():
            x_offset_mm = roi_data["offset_mm"]
            y = self.open_image.center.y 
            y += roi_data["offset_y_mm"] * self.open_image.dpmm
            x_offset_pixels = x_offset_mm * self.open_image.dpmm
            x = x_field_center + x_offset_pixels
            points.append(Point(x, y))
        return points


        

class CustomCTP404(CTP404CP504):
    roi_dist_mm = 58.7  # Default value
    roi_radius_mm = 4  # Smaller diameter
    roi_settings = {
        'Air': {
            'value': AIR,
            'angle': -90,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'PMP': {
            'value': PMP,
            'angle': -120,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'LDPE': {
            'value': LDPE,
            'angle': 180,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'Poly': {
            'value': POLY,
            'angle': 120,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'Acrylic': {
            'value': ACRYLIC,
            'angle': 60,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'Delrin': {
            'value': DELRIN,
            'angle': 0,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
        'Teflon': {
            'value': TEFLON,
            'angle': -60,
            'distance': roi_dist_mm,
            'radius': roi_radius_mm,
        },
    }
    background_roi_settings = {
        '1': {'angle': -30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '2': {'angle': -150, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '3': {'angle': -210, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
        '4': {'angle': 30, 'distance': roi_dist_mm, 'radius': roi_radius_mm},
    }


class CustomCTP486(CTP486):
    '''
    Custom analysis ROIs for the Noise power spectrum analysis
    '''
    nps_rois: dict[str, RectangleROI]
    # Angles are in the polar coordinate system
    nps_roi_settings = {angle: {'value': 0, 'angle': angle, 'distance': 40, 'radius': 15} for angle in range(-90, 260, 20)}
    
    def preprocess(self, catphan) -> None:
        # Save adjacent slice for subtraction in the subsequent nps analysis
        self.adjacent_image = Slice(
            catphan,
            combine_method="mean",
            num_slices=0,
            slice_num=self.slice_num + 1,
            clear_borders=self.clear_borders,
        ).image
    
    def _setup_rois(self) -> None:
        '''Generate our NPS ROIs. They are just square versions of the existing ROIs.'''
        super()._setup_rois()
        self.nps_rois = {}
        for name, setting in self.nps_roi_settings.items():
            self.nps_rois[name] = RectangleROI(
                array=self.image,
                width=setting['radius_pixels'] * 2,
                height=setting['radius_pixels'] * 2,
                angle=setting['angle_corrected'],
                dist_from_center=setting['distance_pixels'],
                phantom_center=self.phan_center,
            )
            
        # ROIs from the adjacent slice, only the array is changed
        self.nps_rois_2 = {}
        for name, setting in self.nps_roi_settings.items():
            self.nps_rois_2[name] = RectangleROI(
                array=self.adjacent_image,
                width=setting['radius_pixels'] * 2,
                height=setting['radius_pixels'] * 2,
                angle=setting['angle_corrected'],
                dist_from_center=setting['distance_pixels'],
                phantom_center=self.phan_center,
            )
            
    @cached_property
    def power_spectrum_2d(self) -> np.ndarray:
        """The power spectrum of the uniformity ROI."""
        return noise_power_spectrum_2d(
            pixel_size=self.mm_per_pixel,
            # Subtract roi from next slice to get only a noise image
            rois=[r.pixel_array - self.nps_rois_2[key].pixel_array 
                  for key, r in self.nps_rois.items()],
        )
    
    @cached_property
    def power_spectrum_1d(self) -> np.ndarray:
        """The 1D power spectrum of the uniformity ROI."""
        return noise_power_spectrum_1d(self.power_spectrum_2d)


# Pass custom modules to the CatPhan model
class CustomCP504(CatPhan504):
    modules = {
        CustomCTP404: {'offset': 0},
        #CTP486: {'offset': -65},
        CustomCTP486: {'offset': -65},
        CTP528CP504: {'offset': 30},
        CTP515: {'offset': -30},
    }
    
# CIRS electron density phantom properties

CIRS_ELECTRON_DENSITY = {
    'Air': {'density': 1.013e-20},
    'Lung(inhale)_IN': {'density': 0.668},
    'Lung(inhale)_OUT': {'density': 0.668},
    'Water': {'density': 1.0},
    'Lung(exhale)_IN': {'density': 1.658},
    'Lung(exhale)_OUT': {'density': 1.658},
    'Adipose_IN': {'density': 3.171},
    'Breast_50/50_IN': {'density': 3.261},
    'Breast_50/50_OUT': {'density': 3.261},
    'Muscle_IN': {'density': 3.483},
    'Muscle_OUT': {'density': 3.483},
    'Liver_IN': {'density': 3.516},
    'Liver_OUT': {'density': 3.516},
    'Bone_200mg/cc_IN': {'density': 3.730},
    'Bone_200mg/cc_OUT': {'density': 3.730},
    'Bone_800mg/cc_IN': {'density': 4.862},
    'Bone_800mg/cc_OUT': {'density': 4.862},
    'Bone_1250mg/cc_OUT': {'density': 5.663},
    'Aluminum_core_IN': {'density': 8.008},
    'Titanium_core_IN': {'density': 12.475},
    'Stainless_steel_IN': {'density': 23.101},
}
CIRS_MASS_DENSITY = {
    'Air': {'density': 0.001204},
    'Lung(inhale)_IN': {'density': 0.205},
    'Lung(inhale)_OUT': {'density': 0.205},
    'Water': {'density': 0.997},
    'Lung(exhale)_IN': {'density': 0.507},
    'Lung(exhale)_OUT': {'density': 0.507},
    'Adipose_IN': {'density': 0.960},
    'Breast_50/50_IN': {'density': 0.990},
    'Breast_50/50_OUT': {'density': 0.990},
    'Muscle_IN': {'density': 1.060},
    'Muscle_OUT': {'density': 1.060},
    'Liver_IN': {'density': 1.070},
    'Liver_OUT': {'density': 1.070},
    'Bone_200mg/cc_IN': {'density': 1.160},
    'Bone_200mg/cc_OUT': {'density': 1.160},
    'Bone_800mg/cc_IN': {'density': 1.530},
    'Bone_800mg/cc_OUT': {'density': 1.530},
    'Bone_1250mg/cc_OUT': {'density': 1.820},
    'Aluminum_core_IN': {'density': 2.700},
    'Titanium_core_IN': {'density': 4.510},
    'Stainless_steel_IN': {'density': 8.030},
}


# Normal density targets

    
class CIRSHUModule(CheeseModule):
    """The pluggable module with user-accessible holes.

    The ROIs of each circle are ~45 degrees apart.
    """

    common_name = "CIRS electron density"
    outer_radius_mm = 115
    outside_radius_mm = 200
    inner_radius_mm = 60
    roi_radius_mm = 9
    metal_radius_mm = 2
    roi_settings = {
        'Air': {
            'angle': -45,
            'distance': outside_radius_mm,
            'radius': roi_radius_mm,
        },
        'Aluminum_core_IN': {
            'angle': 0,
            'distance': 0,
            'radius': metal_radius_mm,
        },
        'Water': {
            'angle': 0,
            'distance': 0,
            'radius': roi_radius_mm - 2,
        },
        'Breast_50/50_IN': {
            'angle': -90,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Muscle_OUT': {
            'angle': -90,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Liver_IN': {
            'angle': -45,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_800mg/cc_OUT': {
            'angle': -45,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(exhale)_IN': {
            'angle': 0,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(inhale)_OUT': {
            'angle': 0,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Adipose_IN': {
            'angle': 45,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Breast_50/50_OUT': {
            'angle': 45,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_800mg/cc_IN': {
            'angle': 90,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        # this one is closer to the ring; presumably because the bottom of the phantom is flatter than the top
        'Liver_OUT': {
            'angle': 90,
            'distance': outer_radius_mm - 5,
            'radius': roi_radius_mm,
        },
        'Muscle_IN': {
            'angle': 135,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_200mg/cc_OUT': {
            'angle': 135,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(inhale)_IN': {
            'angle': 180,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(exhale)_OUT': {
            'angle': 180,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_200mg/cc_IN': {
            'angle': -135,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_1250mg/cc_OUT': {
            'angle': -135,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Titanium_core_IN': {
            'angle': 90,
            'distance': inner_radius_mm,
            'radius': metal_radius_mm,
        },
        'Stainless_steel_IN': {
            'angle': -45,
            'distance': inner_radius_mm,
            'radius': metal_radius_mm,
        },
    }
    
    def _setup_rois(self) -> None:
        # unlike its super, we use simple disk ROIs as we're not doing complicated things.
        for name, setting in self.roi_settings.items():
            self.rois[name] = HighContrastDiskROI(
                self.image,
                setting["angle_corrected"],
                setting["radius_pixels"],
                setting["distance_pixels"],
                self.phan_center,
                contrast_threshold = 0.1
            )
    
    
class EDCIRSModule(CIRS062M):
    module_class = CIRSHUModule
    module: CIRSHUModule
    
    def save_analyzed_image(self, filename: str | Path | BinaryIO, **kwargs) -> None:
        """Save the analyzed summary plot.

        Parameters
        ----------
        filename : str, file object
            The name of the file to save the image to.
        kwargs :
            Any valid matplotlib kwargs.
        """
        self.plot_analyzed_image(show=False, figsize=(10, 10))
        plt.savefig(filename, **kwargs)
        
    def _generate_results_data(self) -> CheeseResult:
        # Get DiskROI attributes normally
        mean_rois = {name: roi.as_dict() for name, roi in self.module.rois.items()}
        # Update with mean value of the ROI
        for name, roi in self.module.rois.items():
            mean_rois[name].update({'mean': roi.mean})
        
        return CheeseResult(
            origin_slice=self.origin_slice,
            num_images=self.num_images,
            phantom_roll=self.catphan_roll,
            rois=mean_rois,
        )
        
# Localization for metal ROIs is not working yet

# =============================================================================
#     def find_origin_slice(self) -> int:
#         """Using a brute force search of the images, find the median HU linearity slice.
# 
#         This method walks through all the images and takes a collapsed circle profile where the HU
#         linearity ROIs are. If the profile contains both low (<800) and high (>800) HU values and most values are the same
#         (i.e. it's not an artifact), then
#         it can be assumed it is an HU linearity slice. The median of all applicable slices is the
#         center of the HU slice.
# 
#         Returns
#         -------
#         int
#             The middle slice of the HU linearity module.
#         """
#         hu_slices = []
#         for image_number in range(0, self.num_images, 2):
#             slice = Slice(
#                 self, image_number, combine=False, clear_borders=self.clear_borders
#             )
#             # print(image_number)
#             # slice.image.plot()
#             if slice.is_phantom_in_view():
#                 circle_prof = CollapsedCircleProfile(
#                     slice.phan_center,
#                     radius=self.localization_radius / self.mm_per_pixel,
#                     image_array=slice.image,
#                     width_ratio=0.05,
#                     num_profiles=5,
#                 )
#                 prof = circle_prof.values
#                 # determine if the profile contains both low and high values and that most values are the same
#                 low_end, high_end = np.percentile(prof, [2, 98])
#                 median = np.median(prof)
#                 middle_variation = np.percentile(prof, 80) - np.percentile(prof, 20)
#                 variation_limit = max(
#                     100, self.dicom_stack.metadata.SliceThickness * -100 + 300
#                 )
#                 if (
#                     (low_end < median - self.hu_origin_slice_variance)
#                     and (high_end > median + self.hu_origin_slice_variance)
#                     and (middle_variation < variation_limit)
#                 ):
#                     hu_slices.append(image_number)
#         
#     def find_phantom_axis(self) -> (Callable, Callable):
#         """We fit all the center locations of the phantom across all slices to a 1D poly function instead of finding them individually for robustness.
# 
#         Normally, each slice would be evaluated individually, but the RadMachine jig gets in the way of
#         detecting the HU module (🤦‍♂️). To work around that in a backwards-compatible way we instead
#         look at all the slices and if the phantom was detected, capture the phantom center.
#         ALL the centers are then fitted to a 1D poly function and passed to the individual slices.
#         This way, even if one slice is messed up (such as because of the phantom jig), the poly function
#         is robust to give the real center based on all the other properly-located positions on the other slices.
#         """
#         z = []
#         center_x = []
#         center_y = []
#         for idx, img in enumerate(self.dicom_stack):
#             slice = Slice(
#                 self,
#                 slice_num=idx,
#                 clear_borders=self.clear_borders,
#                 original_image=img,
#             )
#             if slice.is_phantom_in_view():
#                 roi = slice.phantom_roi
#                 z.append(idx)
#                 center_y.append(roi.centroid[0])
#                 center_x.append(roi.centroid[1])
#         # clip to exclude any crazy values
#         zs = np.array(z)
#         center_xs = np.array(center_x)
#         center_ys = np.array(center_y)
#         # gives an absolute and relative range so tight ranges are all included
#         # but extreme values are excluded. Sometimes the range is very tight
#         # and thus percentiles are not a sure thing
#         x_idxs = np.argwhere(
#             np.isclose(np.median(center_xs), center_xs, atol=3, rtol=0.01)
#         )
#         y_idxs = np.argwhere(
#             np.isclose(np.median(center_ys), center_ys, atol=3, rtol=0.01)
#         )
#         common_idxs = np.intersect1d(x_idxs, y_idxs)
#         # fit to 1D polynomials; inspiration: https://stackoverflow.com/a/45351484
#         fit_zx = np.poly1d(np.polyfit(zs[common_idxs], center_xs[common_idxs], deg=1))
#         fit_zy = np.poly1d(np.polyfit(zs[common_idxs], center_ys[common_idxs], deg=1))
#         return fit_zx, fit_zy
# =============================================================================

    

# Including metal targets

class CIRSMetalHUModule(CheeseModule):
    """The pluggable module with user-accessible holes.

    The ROIs of each circle are ~45 degrees apart.
    """

    common_name = "CIRS electron density, metal ROIs"
    outer_radius_mm = 115
    inner_radius_mm = 60
    roi_radius_mm = 9
    metal_radius_mm = 2
    roi_settings = {
        'Aluminum_core_IN': {
            'angle': 0,
            'distance': 0,
            'radius': metal_radius_mm,
        },
        'Breast_50/50_IN': {
            'angle': -90,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Muscle_OUT': {
            'angle': -90,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Stainless_steel_IN': {
            'angle': -45,
            'distance': inner_radius_mm,
            'radius': metal_radius_mm,
        },
        'Bone_800mg/cc_OUT': {
            'angle': -45,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(exhale)_IN': {
            'angle': 0,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(inhale)_OUT': {
            'angle': 0,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Adipose_IN': {
            'angle': 45,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Breast_50/50_OUT': {
            'angle': 45,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Titanium_core_IN': {
            'angle': 90,
            'distance': inner_radius_mm,
            'radius': metal_radius_mm,
        },
        # this one is closer to the ring; presumably because the bottom of the phantom is flatter than the top
        'Liver_OUT': {
            'angle': 90,
            'distance': outer_radius_mm - 5,
            'radius': roi_radius_mm,
        },
        'Muscle_IN': {
            'angle': 135,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_200mg/cc_OUT': {
            'angle': 135,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(inhale)_IN': {
            'angle': 180,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Lung(exhale)_OUT': {
            'angle': 180,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_200mg/cc_IN': {
            'angle': -135,
            'distance': inner_radius_mm,
            'radius': roi_radius_mm,
        },
        'Bone_1250mg/cc_OUT': {
            'angle': -135,
            'distance': outer_radius_mm,
            'radius': roi_radius_mm,
        },
    }
    
class MetalEDCIRSModule(CIRS062M):
    module_class = CIRSHUModule
    module: CIRSHUModule
    
    def save_analyzed_image(self, filename: str | Path | BinaryIO, **kwargs) -> None:
        """Save the analyzed summary plot.

        Parameters
        ----------
        filename : str, file object
            The name of the file to save the image to.
        kwargs :
            Any valid matplotlib kwargs.
        """
        self.plot_analyzed_image(show=False, figsize=(10, 10))
        plt.savefig(filename, **kwargs)
        
class Slice:
    """Base class for analyzing specific slices of a CBCT dicom set."""

    def __init__(
        self,
        catphan,
        slice_num: int | None = None,
        combine: bool = True,
        combine_method: str = "mean",
        num_slices: int = 0,
        clear_borders: bool = True,
        original_image: ImageLike | None = None,
    ):
        """
        Parameters
        ----------

        catphan : :class:`~pylinac.cbct.CatPhanBase` instance.
            The catphan instance.
        slice_num : int
            The slice number of the DICOM array desired. If None, will use the ``slice_num`` property of subclass.
        combine : bool
            If True, combines the slices +/- ``num_slices`` around the slice of interest to improve signal/noise.
        combine_method : {'mean', 'max'}
            How to combine the slices if ``combine`` is True.
        num_slices : int
            The number of slices on either side of the nominal slice to combine to improve signal/noise; only
            applicable if ``combine`` is True.
        clear_borders : bool
            If True, clears the borders of the image to remove any ROIs that may be present.
        original_image : :class:`~pylinac.core.image.Image` or None
            The array of the slice. This is a bolt-on parameter for optimization.
            Leaving as None is fine, but can increase analysis speed if 1) this image is passed and
            2) there is no combination of slices happening, which is most of the time.
        """
        if slice_num is not None:
            self.slice_num = slice_num
        if combine and num_slices > 0:
            array = combine_surrounding_slices(
                catphan.dicom_stack,
                self.slice_num,
                mode=combine_method,
                slices_plusminus=num_slices,
            )
        elif original_image is not None:
            array = original_image
        else:
            array = catphan.dicom_stack[self.slice_num].array
        self.image = image.load(array)
        self.catphan_size = catphan.catphan_size
        self.mm_per_pixel = catphan.mm_per_pixel
        self.clear_borders = clear_borders
        if catphan._phantom_center_func:
            self._phantom_center_func = catphan._phantom_center_func

    @property
    def __getitem__(self, item):
        return self.image.array[item]

    @cached_property
    def phantom_roi(self) -> RegionProperties:
        """Get the Scikit-Image ROI of the phantom

        The image is analyzed to see if:
        1) the CatPhan is even in the image (if there were any ROIs detected)
        2) an ROI is within the size criteria of the catphan
        3) the ROI area that is filled compared to the bounding box area is close to that of a circle
        """
        # convert the slice to binary and label ROIs
        edges = filters.scharr(self.image.as_type(float))
        if np.max(edges) < 0.1:
            raise ValueError(
                "No edges were found in the image that look like the phantom"
            )
        larr, regionprops, num_roi = get_regions(
            self, fill_holes=True, threshold="mean", clear_borders=self.clear_borders
        )
        # check that there is at least 1 ROI
        if num_roi < 1 or num_roi is None:
            raise ValueError(
                f"The number of ROIs detected {num_roi} was not the number expected (1)"
            )
        catphan_region = sorted(
            regionprops, key=lambda x: np.abs(x.filled_area - self.catphan_size)
        )[0]
        if (self.catphan_size * 1.3 < catphan_region.filled_area) or (
            catphan_region.filled_area < self.catphan_size / 1.3
        ):
            raise ValueError("Unable to find ROI of expected size of the phantom")
        return catphan_region

    def is_phantom_in_view(self) -> bool:
        """Whether the phantom appears to be within the slice."""
        try:
            self.phantom_roi
            return True
        except ValueError:
            return False

    @property
    def phan_center(self) -> Point:
        """Determine the location of the center of the phantom."""
        x = self._phantom_center_func[0](self.slice_num)
        y = self._phantom_center_func[1](self.slice_num)
        return Point(x=x, y=y)
    
def get_regions(
    slice_or_arr: Slice | np.array,
    fill_holes: bool = False,
    clear_borders: bool = True,
    threshold: str = "otsu",
) -> tuple[np.array, list, int]:
    """Get the skimage regions of a black & white image."""
    if threshold == "otsu":
        thresmeth = filters.threshold_otsu
    elif threshold == "mean":
        thresmeth = np.mean
    if isinstance(slice_or_arr, Slice):
        edges = filters.scharr(slice_or_arr.image.array.astype(float))
        center = slice_or_arr.image.center
    elif isinstance(slice_or_arr, np.ndarray):
        edges = filters.scharr(slice_or_arr.astype(float))
        center = (int(edges.shape[1] / 2), int(edges.shape[0] / 2))
    edges = filters.gaussian(edges, sigma=1)
    if isinstance(slice_or_arr, Slice):
        radius = 110 / slice_or_arr.mm_per_pixel
        rr, cc = draw.disk(
            center=(center.y, center.x), radius=radius, shape=edges.shape
        )
        thres = thresmeth(edges[rr, cc])
    else:
        thres = thresmeth(edges)
    bw = edges > thres
    if clear_borders:
        bw = segmentation.clear_border(bw, buffer_size=int(max(bw.shape) / 50))
    if fill_holes:
        bw = ndimage.binary_fill_holes(bw)
    labeled_arr, num_roi = measure.label(bw, return_num=True)
    regionprops = measure.regionprops(labeled_arr, edges)
    return labeled_arr, regionprops, num_roi