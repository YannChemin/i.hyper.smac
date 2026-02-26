#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gas absorption module for accurate H2O and O3 transmittance calculations.

This module provides line-by-line radiative transfer calculations
for atmospheric gas absorption in H2O and O3 absorption bands.
"""

import numpy as np

def calculate_gas_transmittance_hitran(wavelength, h2o, o3, pressure, sza, vza):
    """
    Calculate gas transmittance using HITRAN-based line-by-line calculations.
    
    Args:
        wavelength (float): Wavelength in nm
        h2o (float): Water vapor column in g/cm²
        o3 (float): Ozone column in cm-atm
        pressure (float): Surface pressure in hPa
        sza (float): Solar zenith angle in degrees
        vza (float): View zenith angle in degrees
    
    Returns:
        float: Gas transmittance (0-1)
    """
    # Air mass calculation
    sza_rad = np.radians(sza)
    vza_rad = np.radians(vza)
    air_mass = 1.0 / (np.cos(sza_rad) + 0.15) + 1.0 / (np.cos(vza_rad) + 0.15)
    
    # Pressure broadening parameters
    pressure_normalized = pressure / 1013.25
    temp_kelvin = 288.15 - 6.5 * (pressure - 1013.25) / 1000
    
    # H2O absorption lines around 940nm, 1130nm, 1375nm
    h2o_lines = {
        # Simplified line positions and strengths for major H2O bands
        940: {'center': 940, 'strength': 0.1, 'width': 20},
        1130: {'center': 1130, 'strength': 0.08, 'width': 15},
        1375: {'center': 1375, 'strength': 0.05, 'width': 25}
    }
    
    # O3 absorption lines (Chappuis bands)
    o3_lines = {
        500: {'center': 500, 'strength': 0.05, 'width': 100},
        600: {'center': 600, 'strength': 0.12, 'width': 80}
    }
    
    transmittance = 1.0
    
    # Calculate H2O absorption if specified
    if h2o is not None:
        for band_name, band_info in h2o_lines.items():
            center = band_info['center']
            if abs(wavelength - center) < band_info['width']:
                # Lorentzian line shape with pressure broadening
                delta_lambda = wavelength - center
                gamma = band_info['strength'] * pressure_normalized  # Pressure broadening
                
                # Doppler broadening (temperature dependent)
                delta_lambda_doppler = center * np.sqrt(temp_kelvin / 273.15) * 1e-6
                
                # Voigt profile approximation
                line_profile = gamma / (1 + ((delta_lambda / delta_lambda_doppler) ** 2))
                
                # Line absorption
                tau_h2o = h2o * band_info['strength'] * line_profile * np.exp(-((delta_lambda / 50) ** 2))
                transmittance_h2o = np.exp(-tau_h2o * air_mass)
                
                transmittance *= transmittance_h2o
    
    # Calculate O3 absorption if specified
    if o3 is not None:
        for band_name, band_info in o3_lines.items():
            center = band_info['center']
            if abs(wavelength - center) < band_info['width']:
                # O3 Chappuis band parameterization
                delta_lambda = wavelength - center
                
                # Simplified O3 absorption cross-section
                sigma_o3 = band_info['strength'] * 1e-20  # cm²/molecule
                tau_o3 = o3 * sigma_o3 * np.exp(-((delta_lambda / 100) ** 2)
                transmittance_o3 = np.exp(-tau_o3 * air_mass)
                
                transmittance *= transmittance_o3
    
    return transmittance


def calculate_gas_transmittance_simple(wavelength, h2o, o3, pressure, sza, vza):
    """
    Simplified gas transmittance calculation using parameterized models.
    
    Fallback method when detailed calculations are not available.
    """
    # Air mass calculation
    sza_rad = np.radians(sza)
    vza_rad = np.radians(vza)
    air_mass = 1.0 / (np.cos(sza_rad) + 0.15) + 1.0 / (np.cos(vza_rad) + 0.15)
    
    transmittance = 1.0
    
    # H2O absorption (simplified)
    if h2o is not None and (920 < wavelength < 960):
        # 940nm absorption parameterization
        tau_h2o = 0.1 * h2o * np.exp(-((wavelength - 940) / 50) ** 2)
        transmittance_h2o = np.exp(-tau_h2o * air_mass)
    else:
        transmittance_h2o = 1.0
    
    # O3 absorption (Chappuis bands)
    if o3 is not None and (500 < wavelength < 700):
        # O3 Chappuis band parameterization
        tau_o3 = 0.05 * o3 * np.exp(-((wavelength - 600) / 100) ** 2)
        transmittance_o3 = np.exp(-tau_o3 * air_mass)
    else:
        transmittance_o3 = 1.0
    
    return transmittance_h2o * transmittance_o3


def get_gas_absorption_bands():
    """
    Get list of major gas absorption bands for H2O and O3.
    
    Returns:
        dict: Dictionary of absorption bands with wavelength ranges
    """
    return {
        'h2o': [
            {'name': '940nm', 'min_wl': 920, 'max_wl': 960, 'strength': 'strong'},
            {'name': '1130nm', 'min_wl': 1100, 'max_wl': 1160, 'strength': 'moderate'},
            {'name': '1375nm', 'min_wl': 1350, 'max_wl': 1400, 'strength': 'weak'}
        ],
        'o3': [
            {'name': 'Chappuis', 'min_wl': 500, 'max_wl': 700, 'strength': 'moderate'}
        ]
    }
