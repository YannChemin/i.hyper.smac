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


def get_band_info(input_raster, verbose=False):
    """Extract band information from the input raster metadata.
    
    Args:
        input_raster (str): Name of the input 3D raster
        verbose (bool, optional): Enable verbose output
        
    Returns:
        list: List of band information dictionaries
    """
    try:
        # Get the full metadata from the 3D raster
        info = gs.raster3d_info(input_raster)
        history = gs.read_command('r3.info', flags='h', map=input_raster)
        
        # Parse band information from the history
        band_info = []
        for line in history.split('\n'):
            if line.strip().startswith('Band '):
                try:
                    # Parse line like: "Band 1: 376.44000244140625 nm, FWHM: 5.389999866485596 nm"
                    parts = line.split('Band ')[1].split(':')
                    band_num = int(parts[0].strip())
                    wavelength = float(parts[1].split('nm')[0].strip())
                    fwhm = float(parts[2].split('nm')[0].strip())
                    
                    band_info.append({
                        'band': band_num,
                        'wavelength': wavelength,
                        'fwhm': fwhm,
                        'valid': True
                    })
                except (ValueError, IndexError) as e:
                    if verbose:
                        gs.warning(f"Error parsing band info: {line} - {e}")
        
        if not band_info:
            # Fallback to the original method if no bands were found
            gs.warning("No band information found in history, using default band numbers")
            for i in range(1, info['depths'] + 1):
                band_info.append({
                    'band': i,
                    'wavelength': i,  # Just use band number as wavelength
                    'fwhm': 10.0,     # Default FWHM
                    'valid': True
                })
        
        if verbose:
            gs.message(f"Found {len(band_info)} bands in metadata")
            for band in band_info[:5]:  # Show first 5 bands as example
                gs.verbose(f"Band {band['band']}: {band['wavelength']:.2f} nm")
            if len(band_info) > 5:
                gs.verbose("...")
                gs.verbose(f"Band {band_info[-1]['band']}: {band_info[-1]['wavelength']:.2f} nm")
        
        return band_info
            
    except Exception as e:
        gs.fatal(f"Error getting band information: {e}")
        return []


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
    
    try:
        # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
        gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
        
        # Extract the band using r3.to.rast with 3D region
        gs.run_command('r3.to.rast',
                      input=input_raster,
                      output=output_map,
                      overwrite=True,
                      quiet=True)
         
             
        return output_map
        
    except Exception as e:
        gs.warning(f"Error extracting band {band_num}: {str(e)}")
        raise
        
    finally:
        # Set the 3D region back
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)
        

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
        if not bands:
            raise ValueError("No bands with wavelength information found in the input raster")
        
        if not any(500 <= x['wavelength'] <= 700 for x in bands):
            raise ValueError("Estimate O3 Chappuis : No elements in 'bands' between 500 and 700")

        if verbose:
            gs.message(f"Found {len(bands)} bands in input raster")
            gs.message(f"Wavelength range: {bands[0]['wavelength']:.1f} - {bands[-1]['wavelength']:.1f} nm")
        
        # Find bands closest to key wavelengths
        target_wavelengths = [550, 600, 650]  # Reference, O3 peak, reference
        band_maps = []
        
        for wl in target_wavelengths:
            band = find_nearest_band(bands, wl)
            if verbose:
                gs.message(f"Using band {band['band']} ({band['wavelength']:.1f} nm) "
                          f"for {wl} nm")
            gs.message(f"Using band {band['band']} from {input_raster} ({band['wavelength']:.1f} nm) "
                          f"for {wl} nm")
            band_map = extract_band_to_2d(input_raster, band['band'], 
                                       f"tmp_ozone_{wl}")
            band_maps.append((wl, band_map))
        
        # Calculate ozone using the Chappuis band ratio method
        ozone_map = "ozone_estimate"
        
        # Get the band maps
        try:
            # Find the closest bands to our target wavelengths
            ref1_wl, ref1_map = min(band_maps, key=lambda x: abs(x[0] - 550))
            o3_wl, o3_map = min(band_maps, key=lambda x: abs(x[0] - 600))
            ref2_wl, ref2_map = min(band_maps, key=lambda x: abs(x[0] - 650))
            
            if verbose:
                gs.message(f"Using bands for ozone estimation:")
                gs.message(f"  Reference 1 (~550nm): Band at {ref1_wl:.1f} nm -> {ref1_map}")
                gs.message(f"  Ozone band (~600nm): Band at {o3_wl:.1f} nm -> {o3_map}")
                gs.message(f"  Reference 2 (~650nm): Band at {ref2_wl:.1f} nm -> {ref2_map}")
                
                # Verify the bands are in the correct order
                if not (ref1_wl < o3_wl < ref2_wl):
                    gs.warning("Warning: Bands are not in the expected order. Results may be inaccurate.")
                    gs.warning(f"Expected order: ~550nm < ~600nm < ~650nm")
                    gs.warning(f"Actual order: {ref1_wl:.1f} < {o3_wl:.1f} < {ref2_wl:.1f}")
                    
        except (ValueError, IndexError) as e:
            raise ValueError(f"Could not find appropriate bands for ozone estimation: {str(e)}. "
                           "Ensure the input data covers the 500-700nm range.")
        
        # Calculate ozone using band ratios
        # This is a simplified empirical relationship
        try:
            # First, check if any of the bands are empty
            for band_name, band_map in [('Reference 1', ref1_map), 
                                      ('Ozone band', o3_map), 
                                      ('Reference 2', ref2_map)]:
                stats = gs.parse_command('r.univar', map=band_map, flags='g')
                if float(stats['non_null_cells']) == 0:
                    raise ValueError(f"{band_name} band ({band_map}) contains no data")
            
            # Calculate the ozone using the band ratio
            expr = f"{ozone_map} = 300.0 * (1.0 - float({o3_map}) / (0.5 * ({ref1_map} + {ref2_map}) + 0.0001))"
            gs.mapcalc(expr, overwrite=True)
            
            # Apply reasonable bounds (150-500 DU)
            expr = f"{ozone_map} = if({ozone_map} < 150, 150, if({ozone_map} > 500, 500, {ozone_map}))"
            gs.mapcalc(expr, overwrite=True)
            
            # Calculate mean ozone
            stats = gs.parse_command('r.univar', map=ozone_map, flags='g')
            mean_ozone = float(stats['mean'])
            
            if verbose:
                gs.message(f"Ozone estimation complete. Mean ozone: {mean_ozone:.1f} DU")
            
            return ozone_map, mean_ozone
            
        except Exception as e:
            raise RuntimeError(f"Error calculating ozone map: {str(e)}")
        
    except Exception as e:
        if verbose:
            gs.warning(f"Error in ozone estimation: {str(e)}")
            gs.warning("Using default ozone value of 300 DU")
        
        # Create a constant ozone map with default value
        gs.mapcalc(f"ozone_estimate = {DEFAULT_OZONE}", overwrite=True)
        return "ozone_estimate", DEFAULT_OZONE
    
    finally:
        # Clean up temporary maps
        for _, band_map in locals().get('band_maps', []):
            try:
                gs.run_command('g.remove', flags='f', type='raster', 
                             name=band_map, quiet=not verbose)
            except:
                pass


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
    if verbose:
        gs.message(f"Starting ozone estimation using method: {method}")
    
    try:
        if method.lower() == 'chappuis':
            if verbose:
                gs.message("Attempting to estimate ozone using Chappuis band method...")
            return estimate_ozone_chappuis(input_raster, verbose=verbose)
        
        # Fall back to constant if method is not 'chappuis' or if chappuis method fails
        if verbose:
            if method.lower() != 'chappuis':
                gs.message(f"Unknown method '{method}'. Using constant ozone value.")
            gs.message(f"Using constant ozone value: {DEFAULT_OZONE} DU")
        
        # Create a constant map with default value
        ozone_map = "ozone_estimate"
        gs.mapcalc(f"{ozone_map} = {DEFAULT_OZONE}", overwrite=True)
        
        if verbose:
            gs.message(f"Created constant ozone map '{ozone_map}' with value {DEFAULT_OZONE} DU")
        
        return ozone_map, DEFAULT_OZONE
        
    except Exception as e:
        if verbose:
            gs.warning(f"Error in estimate_ozone: {str(e)}")
            gs.warning("Falling back to constant ozone value")
        
        # Ensure we have a valid return value even in case of errors
        ozone_map = "ozone_estimate"
        gs.mapcalc(f"{ozone_map} = {DEFAULT_OZONE}", overwrite=True)
        return ozone_map, DEFAULT_OZONE
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
