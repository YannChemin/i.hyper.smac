#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Water Vapor Content (WVC) estimation module for hyperspectral imagery.

This module provides functions to estimate Water Vapor Content from
hyperspectral data using atmospheric absorption features, primarily
around 940nm and 1130nm water vapor absorption bands.
"""

import os
import numpy as np
import grass.script as gs

# Water vapor absorption features (in nm)
WATER_VAPOR_BANDS = {
    '940': {
        'center': 940,
        'window': (900, 980),  # Window for continuum removal
        'absorption': (930, 950)  # Absorption feature range
    },
    '1130': {
        'center': 1130,
        'window': (1080, 1180),
        'absorption': (1120, 1140)
    }
}

# Default parameters for WVC estimation
DEFAULT_WVC_MIN = 0.1  # g/cm²
DEFAULT_WVC_MAX = 8.0  # g/cm²


def get_raster3d_info(raster3d):
    """Get information about 3D raster."""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")


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


def find_nearest_band(band_info, target_wavelength):
    """Find the band closest to the target wavelength.
    
    Args:
        band_info (list): List of band information dictionaries
        target_wavelength (float): Target wavelength in nm
        
    Returns:
        dict: Band information for the closest band
    """
    return min(band_info, key=lambda x: abs(x['wavelength'] - target_wavelength))

def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers."""
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


class WVCEstimator:
    """Class for estimating Water Vapor Content from hyperspectral data."""
    
    def __init__(self, input_raster, dem, sensor_config=None, verbose=False):
        """Initialize the WVC estimator.
        
        Args:
            input_raster (str): Name of the input 3D hyperspectral raster
            dem (str): Name of the digital elevation model (DEM) raster
            sensor_config (dict, optional): Sensor configuration dictionary
            verbose (bool, optional): Enable verbose output
        """
        self.input_raster = input_raster
        self.dem = dem
        self.sensor_config = sensor_config or {}
        self.verbose = verbose
        self.band_info = []
        self.wvc_map = None
        self.temp_maps = []  # Track temporary maps for cleanup
        
        # Get band information
        self._get_band_info()
    
    def _get_band_info(self):
        """Extract band information from the 3D raster's metadata."""
        #import re
        
        # Get the full metadata from the 3D raster
        #info = gs.raster3d_info(self.input_raster)
        history = gs.read_command('r3.info', flags='h', map=self.input_raster)
        
        # Parse band information from the history
        self.band_info = []
        for line in history.split('\n'):
            if line.strip().startswith('Band '):
                try:
                    # Parse line like: "Band 1: 376.44000244140625 nm, FWHM: 5.389999866485596 nm"
                    parts = line.split('Band ')[1].split(':')
                    band_num = int(parts[0].strip())
                    wavelength = float(parts[1].split('nm')[0].strip())
                    fwhm = float(parts[2].split('nm')[0].strip())
                    
                    self.band_info.append({
                        'band': band_num,
                        'wavelength': wavelength,
                        'fwhm': fwhm,
                        'valid': True
                    })
                except (ValueError, IndexError) as e:
                    if self.verbose:
                        gs.warning(f"Error parsing band info: {line} - {e}")
        
        if not self.band_info:
            gs.fatal(f"No wavelength metadata found in {self.input_raster}")
        
        if self.verbose:
            gs.message(f"Found {len(self.band_info)} bands in metadata")
            for band in self.band_info[:5]:  # Show first 5 bands as example
                gs.verbose(f"Band {band['band']}: {band['wavelength']:.2f} nm")
            if len(self.band_info) > 5:
                gs.verbose("...")
                gs.verbose(f"Band {self.band_info[-1]['band']}: {self.band_info[-1]['wavelength']:.2f} nm")
            
    def _find_bands_in_range(self, min_wl, max_wl):
        """Find bands within a wavelength range.
        
        Args:
            min_wl (float): Minimum wavelength (nm)
            max_wl (float): Maximum wavelength (nm)
            
        Returns:
            list: List of band dictionaries within the range
        """
        return [b for b in self.band_info if min_wl <= b['wavelength'] <= max_wl]
    
    def _get_band_by_wavelength(self, target_wl):
        """Get the band closest to the target wavelength.
        
        Args:
            target_wl (float): Target wavelength in nm
            
        Returns:
            dict: Band information for the closest band
        """
        return min(self.band_info, key=lambda x: abs(x['wavelength'] - target_wl))

    def parse_wavelength_from_metadata(self, band_num):
        """Parse wavelength from band metadata by extracting to temporary 2D raster."""
        temp_band = f"tmp_wvc_meta_{os.getpid()}_{band_num}"
        
        try:
            temp_band = self._extract_band_to_2d(
                band_num=band_num,
                output_map=temp_band,
            )
            
            # Read metadata
            result = gs.read_command('r.support', map=temp_band, flags='n')
            
            wavelength = None
            fwhm = None
            valid = True
            unit = "nm"
            
            for line in result.split('\n'):
                line = line.strip()
                if line.startswith('wavelength='):
                    wavelength = float(line.split('=')[1])
                elif line.startswith('FWHM='):
                    fwhm = float(line.split('=')[1])
                elif line.startswith('valid='):
                    valid = int(line.split('=')[1]) == 1
                elif line.startswith('unit='):
                    unit = line.split('=')[1].strip()
            
            return wavelength, fwhm, valid, unit
            
        except Exception as e:
            gs.warning(f"Could not read metadata for band {band_num}: {e}")
            return None, None, True, "nm"
    
        finally:
            # Clean up temporary band
            gs.run_command('g.remove', flags='f', type='raster', name=temp_band, quiet=True)
            # Also clean up any temporary maps that might have been created
            gs.run_command('g.remove', flags='f', type='raster_3d', pattern=f'tmp_*', quiet=True)
            gs.run_command('g.remove', flags='f', type='raster', pattern=f'tmp_band_*', quiet=True)


    def _extract_band_to_2d(self, band_num, output_map=None):
        """Extract a single band from 3D raster to 2D raster using g.region and r3.to.rast.
        
        Args:
            band_num (int): Band number to extract (1-based index)
            output_map (str, optional): Name for the output 2D raster. 
                                      If None, a temporary name is generated.
            
        Returns:
            str: Name of the extracted 2D raster
        """
        import time
        import random
        
        if not output_map:
            output_map = f"tmp_band_{band_num}_{int(time.time())}_{random.randint(1000, 9999)}"
        
        try:
            # Clean up any existing files with the same name
            gs.run_command('g.remove', flags='f', type='raster', pattern=f"{output_map}*", quiet=True)
            
            # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
            
            # Convert the 3D raster to 2D with overwrite
            gs.run_command('r3.to.rast',
                         input=self.input_raster,
                         output=output_map,
                         overwrite=True,
                         quiet=not self.verbose)
            
            # The output will be named output_map_00001
            output_file = f"{output_map}_00001"
            
            # Rename the output file to the desired name
            gs.run_command('g.rename',
                         raster=f"{output_file},{output_map}",
                         overwrite=True,
                         quiet=True)
            
            if self.verbose:
                gs.message(f"Extracted band {band_num} to {output_map}")
            
            # Add to temporary maps for cleanup
            self.temp_maps.append(output_map)
            self.temp_maps.append(output_file)
            
            return output_map
            
        except Exception as e:
            gs.warning(f"Error extracting band {band_num}: {e}")
            raise
            
        finally:
            # Clean up any temporary files
            gs.run_command('g.remove', flags='f', type='raster', pattern=f"{output_map}_00001", quiet=True)
            gs.run_command('g.remove', flags='f', type='raster_3d', name='RASTER3D_MASK', quiet=True)
            
    def _calculate_ndvi(self):
        """Calculate NDVI for vegetation masking."""
        # Find nearest bands for red and NIR
        red_band = self._get_band_by_wavelength(670)  # Red band ~670nm
        nir_band = self._get_band_by_wavelength(870)  # NIR band ~870nm
        
        if self.verbose:
            gs.message(f"Calculating NDVI using bands:")
            gs.message(f"  Red: {red_band['wavelength']:.1f} nm (band {red_band['band']})")
            gs.message(f"  NIR: {nir_band['wavelength']:.1f} nm (band {nir_band['band']})")
        
        # Extract bands from 3D raster
        red_map = self._extract_band_to_2d(red_band['band'])
        nir_map = self._extract_band_to_2d(nir_band['band'])
        self.temp_maps.extend([red_map, nir_map])
        
        # Calculate NDVI: (NIR - Red) / (NIR + Red)
        ndvi_map = f"tmp_wvc_ndvi_{os.getpid()}"
        expr = f"{ndvi_map} = float({nir_map} - {red_map}) / float({nir_map} + {red_map} + 0.0001)"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        
        self.temp_maps.append(ndvi_map)
        return ndvi_map
    
    def _estimate_wvc_absorption_feature(self, feature_name):
        """Estimate WVC using a specific water vapor absorption feature.
        
        Args:
            feature_name (str): Name of the absorption feature ('940' or '1130')
            
        Returns:
            tuple: (wvc_map, mean_wvc) where wvc_map is the WVC raster and mean_wvc is the mean value
        """
        feature = WATER_VAPOR_BANDS[feature_name]
        center_wl = feature['center']
        min_wl, max_wl = feature['window']
        abs_min, abs_max = feature['absorption']
        
        # Find bands in the window
        window_bands = self._find_bands_in_range(min_wl, max_wl)
        if len(window_bands) < 3:
            raise ValueError(f"Not enough bands in the {center_wl}nm window "
                           f"(found {len(window_bands)}, need at least 3)")
        
        # Find absorption band (closest to center wavelength)
        abs_band = min(window_bands, key=lambda x: abs(x['wavelength'] - center_wl))
        
        # Find continuum bands (outside absorption feature)
        left_bands = [b for b in window_bands if b['wavelength'] < abs_min]
        right_bands = [b for b in window_bands if b['wavelength'] > abs_max]
        
        if not left_bands or not right_bands:
            raise ValueError(f"Could not find continuum bands for {center_wl}nm feature. "
                           f"Left: {len(left_bands)}, Right: {len(right_bands)}")
        
        # Use the closest bands to the absorption feature for continuum
        left_band = max(left_bands, key=lambda x: x['wavelength'])
        right_band = min(right_bands, key=lambda x: x['wavelength'])
        
        if self.verbose:
            gs.message(f"Using bands for {center_wl}nm WVC estimation:")
            gs.message(f"  Left continuum: {left_band['wavelength']:.1f}nm (band {left_band['band']})")
            gs.message(f"  Absorption: {abs_band['wavelength']:.1f}nm (band {abs_band['band']})")
            gs.message(f"  Right continuum: {right_band['wavelength']:.1f}nm (band {right_band['band']})")
        
        # Extract bands from 3D raster
        gs.message(f"Extracting bands for {center_wl}nm WVC estimation:")
        left_map = self._extract_band_to_2d(left_band['band'])
        gs.message(f"Left band: {left_map}")
        abs_map = self._extract_band_to_2d(abs_band['band'])
        gs.message(f"Absorption band: {abs_map}")
        right_map = self._extract_band_to_2d(right_band['band'])
        gs.message(f"Right band: {right_map}")
        self.temp_maps.extend([left_map, abs_map, right_map])
        
        # Create a mask for valid pixels (e.g., clear land, not clouds or water)
        ndvi_map = self._calculate_ndvi()
        valid_mask = f"tmp_wvc_valid_mask_{os.getpid()}"
        # Use pixels with moderate NDVI (land, not water or clouds)
        expr = f"{valid_mask} = if({ndvi_map} > 0.1 && {ndvi_map} < 0.9, 1, null())"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(valid_mask)
        
        # Calculate continuum-removed reflectance (band depth)
        # Continuum is a linear interpolation between left and right shoulders
        wl_left = left_band['wavelength']
        wl_abs = abs_band['wavelength']
        wl_right = right_band['wavelength']
        
        # Weight for linear interpolation
        weight = (wl_abs - wl_left) / (wl_right - wl_left)
        
        # Calculate continuum and band depth
        band_depth = f"tmp_wvc_depth_{feature_name}_{os.getpid()}"
        expr = (f"{band_depth} = if({valid_mask}, "
                f"1.0 - ({abs_map} / "
                f"({left_map} + ({right_map} - {left_map}) * {weight} + 0.0001)), "
                f"null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(band_depth)
        
        # Convert band depth to WVC (g/cm²)
        # This uses an empirical relationship between band depth and water vapor
        # The relationship depends on the specific absorption feature
        if feature_name == '940':
            # 940nm feature: stronger absorption, more sensitive
            scale_factor = 3.5
            offset = 0.2
        else:  # 1130nm
            # 1130nm feature: weaker absorption, less sensitive
            scale_factor = 5.0
            offset = 0.5
        
        wvc_map = f"tmp_wvc_{feature_name}_{os.getpid()}"
        expr = f"{wvc_map} = if({valid_mask}, {band_depth} * {scale_factor} + {offset}, null())"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(wvc_map)
        
        # Get statistics
        try:
            stats = gs.parse_command('r.univar', map=wvc_map, flags='g')
            mean_wvc = float(stats['mean'])
            
            # Apply min/max constraints
            mean_wvc = max(DEFAULT_WVC_MIN, min(mean_wvc, DEFAULT_WVC_MAX))
            
            if self.verbose:
                gs.message(f"Estimated WVC from {center_wl}nm feature: {mean_wvc:.3f} g/cm²")
            
        except Exception as e:
            gs.warning(f"Failed to calculate WVC statistics: {e}")
            mean_wvc = 2.0  # Default value
        
        return wvc_map, mean_wvc
    
    def estimate_wvc(self, method='average'):
        """Estimate Water Vapor Content.
        
        Args:
            method (str): WVC estimation method ('940nm', '1130nm', or 'average')
            
        Returns:
            tuple: (wvc_map, mean_wvc) where wvc_map is the WVC raster and mean_wvc is the mean value
        """
        if method not in ['940nm', '1130nm', 'average']:
            raise ValueError("method must be '940nm', '1130nm', or 'average'")
        
        wvc_940_map = None
        wvc_1130_map = None
        wvc_940_val = None
        wvc_1130_val = None
        
        # Try 940nm method
        if method in ['940nm', 'average']:
            try:
                wvc_940_map, wvc_940_val = self._estimate_wvc_absorption_feature('940')
                if method == '940nm':
                    self.wvc_map = wvc_940_map
                    return wvc_940_map, wvc_940_val
            except Exception as e:
                if self.verbose:
                    gs.warning(f"940nm WVC estimation failed: {str(e)}")
                if method == '940nm':
                    raise
        
        # Try 1130nm method
        if method in ['1130nm', 'average']:
            try:
                wvc_1130_map, wvc_1130_val = self._estimate_wvc_absorption_feature('1130')
                if method == '1130nm':
                    self.wvc_map = wvc_1130_map
                    return wvc_1130_map, wvc_1130_val
            except Exception as e:
                if self.verbose:
                    gs.warning(f"1130nm WVC estimation failed: {str(e)}")
                if method == '1130nm':
                    raise
        
        # Average method
        if method == 'average':
            if wvc_940_map and wvc_1130_map:
                # Average the two estimates
                wvc_map = f"tmp_wvc_avg_{os.getpid()}"
                expr = f"{wvc_map} = ({wvc_940_map} + {wvc_1130_map}) / 2.0"
                gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
                self.temp_maps.append(wvc_map)
                
                mean_wvc = (wvc_940_val + wvc_1130_val) / 2.0
                self.wvc_map = wvc_map
                
                if self.verbose:
                    gs.message(f"Average WVC from both features: {mean_wvc:.3f} g/cm²")
                
                return wvc_map, mean_wvc
            elif wvc_940_map:
                self.wvc_map = wvc_940_map
                return wvc_940_map, wvc_940_val
            elif wvc_1130_map:
                self.wvc_map = wvc_1130_map
                return wvc_1130_map, wvc_1130_val
        
        # If we get here, all methods failed
        raise RuntimeError("Could not estimate WVC using any available method. "
                         "Check that your data includes NIR bands with water vapor features.")
    
    def cleanup(self):
        """Clean up temporary maps."""
        if self.verbose:
            gs.message("Cleaning up temporary maps...")
        
        for map_name in self.temp_maps:
            if gs.find_file(map_name, element='cell')['file']:
                gs.run_command('g.remove', flags='f', type='raster', 
                             name=map_name, quiet=True)


def estimate_wvc(input_raster, dem, method='average', sensor_config=None, verbose=False):
    """Convenience function to estimate WVC from hyperspectral data.
    
    Args:
        input_raster (str): Name of the input 3D hyperspectral raster
        dem (str): Name of the digital elevation model (DEM) raster
        method (str, optional): WVC estimation method. Defaults to 'average'.
        sensor_config (dict, optional): Sensor configuration dictionary
        verbose (bool, optional): Enable verbose output
        
    Returns:
        tuple: (wvc_map, mean_wvc) where wvc_map is the WVC raster and mean_wvc is the mean value
    """
    estimator = WVCEstimator(input_raster, dem, sensor_config, verbose)
    try:
        wvc_map, mean_wvc = estimator.estimate_wvc(method)
        
        # Create a permanent copy
        permanent_wvc = f"{input_raster}_wvc"
        gs.run_command('g.copy', raster=f"{wvc_map},{permanent_wvc}", 
                      overwrite=True, quiet=not verbose)
        
        return permanent_wvc, mean_wvc
    finally:
        estimator.cleanup()


if __name__ == "__main__":
    # Example usage
    # python wvc.py prisma_hyper dem=dem_10m --method=average --verbose
    import argparse
    
    parser = argparse.ArgumentParser(description='Estimate Water Vapor Content from hyperspectral data')
    parser.add_argument('input', help='Input 3D hyperspectral raster')
    parser.add_argument('dem', help='Digital Elevation Model (DEM) raster')
    parser.add_argument('--method', default='average', 
                      choices=['940nm', '1130nm', 'average'],
                      help='WVC estimation method')
    parser.add_argument('--output', help='Output WVC map name')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    wvc_map, mean_wvc = estimate_wvc(args.input, args.dem, args.method, None, args.verbose)
    
    if args.output and args.output != wvc_map:
        gs.run_command('g.rename', raster=f"{wvc_map},{args.output}", overwrite=True)
        wvc_map = args.output
    
    print(f"WVC map created: {wvc_map}")
    print(f"Mean WVC: {mean_wvc:.3f} g/cm²")
