#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shared utility functions for i.hyper.smac atmospheric correction module.

This module provides common functions used across the library for band
information extraction, wavelength conversion, and raster operations.
"""

import os
import grass.script as gs


def get_band_info(input_raster, verbose=False):
    """Extract band information from the input raster metadata.

    Args:
        input_raster (str): Name of the input 3D raster
        verbose (bool, optional): Enable verbose output

    Returns:
        list: List of band information dictionaries with keys:
            - 'band': Band number (1-based)
            - 'wavelength': Central wavelength in nm
            - 'fwhm': Full width at half maximum in nm
            - 'valid': Whether band data is valid
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
            for i in range(1, int(info['depths']) + 1):
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
    """Find the band with wavelength closest to the target.

    Args:
        bands (list): List of band information dictionaries (from get_band_info)
                     or list of tuples (wavelength, map_name)
        target_wavelength (float): Target wavelength in nm

    Returns:
        dict: Band information for the closest band
    """
    # If bands is a list of tuples (wavelength, map_name), convert to list of dicts
    if bands and isinstance(bands[0], (list, tuple)):
        bands = [{'wavelength': wl, 'band': i + 1} for i, (wl, _) in enumerate(bands)]

    return min(bands, key=lambda x: abs(x['wavelength'] - target_wavelength))


def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers.

    Args:
        wavelength (float): Wavelength value
        unit (str): Unit of the wavelength ('nm', 'um', 'm', etc.)

    Returns:
        float: Wavelength in nanometers
    """
    unit = unit.lower().strip()

    if unit in ['nm', 'nanometer', 'nanometers']:
        return wavelength
    elif unit in ['um', 'µm', 'micrometer', 'micrometers', 'micron', 'microns']:
        return wavelength * 1000.0
    elif unit in ['m', 'meter', 'meters']:
        return wavelength * 1e9
    else:
        gs.warning(f"Unknown wavelength unit '{unit}', assuming nanometers")
        return wavelength


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
        output_map = f"tmp_band_{band_num}_{os.getpid()}"

    try:
        # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
        gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)

        # Extract the band using r3.to.rast with 3D region
        gs.run_command('r3.to.rast',
                       input=input_raster,
                       output=output_map,
                       overwrite=True,
                       quiet=True)

        # The output will be named output_map_00001
        output_file = f"{output_map}_00001"

        # Rename the output file to the desired name
        gs.run_command('g.rename',
                       raster=f"{output_file},{output_map}",
                       overwrite=True,
                       quiet=True)

        return output_map

    except Exception as e:
        gs.warning(f"Error extracting band {band_num}: {str(e)}")
        raise

    finally:
        # Set the 3D region back
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)


def get_raster3d_info(raster3d):
    """Get information about 3D raster.

    Args:
        raster3d (str): Name of the 3D raster

    Returns:
        dict: Dictionary with raster3d information
    """
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")
