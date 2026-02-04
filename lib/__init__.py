"""
Support library for i.hyper.smac atmospheric correction module.

This package provides utilities for atmospheric correction of hyperspectral
imagery using the SMAC (Simplified Method for Atmospheric Correction) algorithm.
"""

# Import shared utility functions
from .utils import (
    get_band_info,
    find_nearest_band,
    convert_wavelength_to_nm,
    extract_band_to_2d,
    get_raster3d_info,
)

# Import submodules for convenience
from .aod import estimate_aod, AODEstimator
from .wvc import estimate_wvc, WVCEstimator
from .o3 import estimate_ozone
from .smac_coef_generator import (
    SMACCoefficients,
    SMACCoefficientGenerator,
    generate_smac_coefficients,
)

__all__ = [
    # Utility functions
    'get_band_info',
    'find_nearest_band',
    'convert_wavelength_to_nm',
    'extract_band_to_2d',
    'get_raster3d_info',
    # Estimation functions
    'estimate_aod',
    'AODEstimator',
    'estimate_wvc',
    'WVCEstimator',
    'estimate_ozone',
    # Coefficient generation
    'SMACCoefficients',
    'SMACCoefficientGenerator',
    'generate_smac_coefficients',
]
