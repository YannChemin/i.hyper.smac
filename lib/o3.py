#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ozone (O₃) estimation module for hyperspectral imagery.

This module provides functions to estimate total column ozone from
hyperspectral data using the Chappuis band absorption feature.
"""

import os
import numpy as np
import grass.script as gs
from pathlib import Path

# Default ozone value (Dobson Units)
DEFAULT_OZONE = 300.0  # Typical mid-latitude value (250-350 DU)

# Ozone absorption coefficients (cm⁻¹/(atm·cm))
# These are approximate values for the Chappuis band
OZONE_ABSORPTION = {
    550: 1.0e-4,  # Reference wavelength
    600: 5.0e-4,  # Ozone absorption peak in Chappuis band
    650: 2.5e-4   # Reference wavelength
}


def get_band_info(raster3d, verbose=False):
    """Extract band information from 3D raster metadata.
    
    Args:
        raster3d (str): Name of the input 3D raster
        verbose (bool, optional): Enable verbose output
        
    Returns:
        list: List of dictionaries with band information
    """
    try:
        # Get band information using r3.info
        info = gs.parse_command('r3.info', flags='g', map=raster3d)
        
        # Extract band metadata (assuming format: Wavelength=xxx.xxnm,FWHM=yy.yynm)
        bands = []
        for i in range(1, int(info['depths']) + 1):
            band_info = {}
            # Get band metadata (wavelength and FWHM)
            metadata = gs.parse_command('r3.info', flags='m', map3d=raster3d, 
                                      depth=str(i))
            
            # Extract wavelength and FWHM from metadata
            for line in metadata.split('\n'):
                if 'Wavelength' in line:
                    band_info['wavelength'] = float(line.split('=')[1].replace('nm', '').strip())
                elif 'FWHM' in line:
                    band_info['fwhm'] = float(line.split('=')[1].replace('nm', '').strip())
            
            if 'wavelength' in band_info:
                band_info['band'] = i
                bands.append(band_info)
                
        if verbose:
            gs.message(f"Found {len(bands)} bands with wavelength information")
            
        return sorted(bands, key=lambda x: x['wavelength'])
        
    except Exception as e:
        gs.fatal(f"Error getting band information: {str(e)}")


def find_nearest_band(bands, target_wavelength):
    """Find the band closest to the target wavelength.
    
    Args:
        bands (list): List of band information dictionaries
        target_wavelength (float): Target wavelength in nm
        
    Returns:
        dict: Band information for the closest band
    """
    return min(bands, key=lambda x: abs(x['wavelength'] - target_wavelength))


def extract_band_to_2d(input_raster, band_num, output_map=None):
    """Extract a single band from 3D raster to 2D raster.
    
    Args:
        input_raster (str): Input 3D raster name
        band_num (int): Band number to extract (1-based index)
        output_map (str, optional): Output 2D raster name
        
    Returns:
        str: Name of the extracted 2D raster
    """
    if not output_map:
        output_map = f"tmp_band_{band_num}"
    
    # Clean up any existing files with the same name
    gs.run_command('g.remove', flags='f', type='raster', 
                  pattern=f"{output_map}*", quiet=True)
    
    # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
    gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
    
    # Extract the band using r3.to.rast
    gs.run_command('r3.to.rast',
                  input=input_raster,
                  output=output_map,
                  overwrite=True,
                  quiet=True)
    
    # Set the 3D region back
    gs.run_command('g.region', raster_3d=input_raster, quiet=True)
    
    return output_map


def estimate_ozone_chappuis(input_raster, verbose=False):
    """Estimate total column ozone using the Chappuis band absorption.
    
    This method uses the differential absorption in the Chappuis band (500-700 nm)
    to estimate total column ozone. It's less accurate than UV methods but works
    with visible-range hyperspectral data.
    
    Args:
        input_raster (str): Input 3D hyperspectral raster
        verbose (bool, optional): Enable verbose output
        
    Returns:
        tuple: (ozone_map, mean_ozone) where ozone_map is the raster name and 
               mean_ozone is the mean value in Dobson Units (DU)
    """
    try:
        # Get band information
        bands = get_band_info(input_raster, verbose=verbose)
        
        # Find bands closest to key wavelengths
        target_wavelengths = [550, 600, 650]  # Reference, O3 peak, reference
        band_maps = []
        
        for wl in target_wavelengths:
            band = find_nearest_band(bands, wl)
            if verbose:
                gs.message(f"Using band {band['band']} ({band['wavelength']:.1f} nm) "
                          f"for {wl} nm")
            band_map = extract_band_to_2d(input_raster, band['band'], 
                                        f"tmp_ozone_{wl}")
            band_maps.append((wl, band_map))
        
        # Calculate ozone using the Chappuis band ratio method
        # This is a simplified approach - for production, consider using a LUT or RTM
        ozone_map = "ozone_estimate"
        
        # Get the band maps
        ref1_band = next(b for wl, b in band_maps if abs(wl - 550) < 10)
        o3_band = next(b for wl, b in band_maps if abs(wl - 600) < 10)
        ref2_band = next(b for wl, b in band_maps if abs(wl - 650) < 10)
        
        # Calculate ozone using band ratios
        # This is a simplified empirical relationship
        expr = f"{ozone_map} = 300.0 * (1.0 - float({o3_band}) / (0.5 * ({ref1_band} + {ref2_band}) + 0.0001))"
        gs.mapcalc(expr, overwrite=True)
        
        # Apply reasonable bounds (150-500 DU)
        expr = f"{ozone_map} = if({ozone_map} < 150, 150, if({ozone_map} > 500, 500, {ozone_map}))"
        gs.mapcalc(expr, overwrite=True)
        
        # Calculate mean ozone
        stats = gs.parse_command('r.univar', map=ozone_map, flags='g')
        mean_ozone = float(stats['mean'])
        
        if verbose:
            gs.message(f"Ozone estimation complete. Mean ozone: {mean_ozone:.1f} DU")
        
        # Clean up temporary maps
        for _, band_map in band_maps:
            gs.run_command('g.remove', flags='f', type='raster', 
                          name=band_map, quiet=True)
        
        return ozone_map, mean_ozone
        
    except Exception as e:
        if verbose:
            gs.warning(f"Error estimating ozone: {str(e)}. Using default value.")
        
        # Create a constant ozone map with default value
        gs.mapcalc(f"ozone_estimate = {DEFAULT_OZONE}", overwrite=True)
        return "ozone_estimate", DEFAULT_OZONE


def estimate_ozone(input_raster, method='chappuis', verbose=False):
    """Estimate total column ozone from hyperspectral data.
    
    Args:
        input_raster (str): Input 3D hyperspectral raster
        method (str): Ozone estimation method ('chappuis' or 'constant')
        verbose (bool, optional): Enable verbose output
        
    Returns:
        tuple: (ozone_map, mean_ozone) where ozone_map is the raster name and 
               mean_ozone is the mean value in Dobson Units (DU)
    """
    if method == 'chappuis':
        return estimate_ozone_chappuis(input_raster, verbose=verbose)
    elif method == 'constant':
        if verbose:
            gs.message(f"Using constant ozone value: {DEFAULT_OZONE} DU")
        gs.mapcalc(f"ozone_estimate = {DEFAULT_OZONE}", overwrite=True)
        return "ozone_estimate", DEFAULT_OZONE
    else:
        raise ValueError(f"Unknown ozone estimation method: {method}")


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) != 2:
        gs.fatal("Usage: python o3.py <input_3d_raster>")
    
    input_raster = sys.argv[1]
    ozone_map, mean_ozone = estimate_ozone(input_raster, verbose=True)
    if gs.verbosity() > 0:
        gs.message(f"Ozone estimation complete. Mean ozone: {mean_ozone:.1f} DU")
        gs.message(f"Ozone map created: {ozone_map}")
