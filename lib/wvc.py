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


def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength from band metadata by extracting to temporary 2D raster."""
    temp_band = f"tmp_wvc_meta_{os.getpid()}_{band_num}"
    
    try:
        # Extract single level from 3D raster
        gs.run_command('r3.to.rast', input=raster3d, output=temp_band,
                      level=band_num, quiet=True, overwrite=True)
        
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
        if gs.find_file(temp_band, element='cell')['file']:
            gs.run_command('g.remove', type='raster', name=temp_band, 
                          flags='f', quiet=True)


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


def extract_band_to_2d(raster3d, band_num, output_name=None):
    """Extract a band from 3D raster to 2D raster.
    
    Args:
        raster3d (str): Name of the 3D raster
        band_num (int): Band number to extract
        output_name (str, optional): Output raster name. If None, generates temp name.
        
    Returns:
        str: Name of the extracted 2D raster
    """
    if output_name is None:
        output_name = f"tmp_wvc_band_{os.getpid()}_{band_num}"
    
    gs.run_command('r3.to.rast', input=raster3d, output=output_name,
                  level=band_num, quiet=True, overwrite=True)
    
    return output_name


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
        """Extract band information from the input raster."""
        info = get_raster3d_info(self.input_raster)
        num_bands = int(info['depths'])
        
        if self.verbose:
            gs.message(f"Scanning {num_bands} bands for wavelength metadata...")
        
        # Get band metadata
        for i in range(1, num_bands + 1):
            wavelength, fwhm, valid, unit = parse_wavelength_from_metadata(
                self.input_raster, i
            )
            
            if wavelength is not None:
                wavelength_nm = convert_wavelength_to_nm(wavelength, unit)
                self.band_info.append({
                    'band': i,
                    'wavelength': wavelength_nm,
                    'fwhm': fwhm if fwhm else 10,
                    'valid': valid
                })
                
                if self.verbose:
                    gs.verbose(f"Band {i}: {wavelength_nm:.2f} nm (valid: {valid})")
        
        if not self.band_info:
            # Fallback: use sensor config or default wavelength range
            gs.warning("No wavelength metadata found, using default wavelength distribution")
            
            if self.sensor_config and 'range' in self.sensor_config:
                min_wl, max_wl = self.sensor_config['range']
            else:
                min_wl, max_wl = 400, 2500
            
            for i in range(1, num_bands + 1):
                wavelength = min_wl + (max_wl - min_wl) * (i - 1) / (num_bands - 1)
                self.band_info.append({
                    'band': i,
                    'wavelength': wavelength,
                    'fwhm': 10,
                    'valid': True
                })
    
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
        red_map = extract_band_to_2d(self.input_raster, red_band['band'])
        nir_map = extract_band_to_2d(self.input_raster, nir_band['band'])
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
        left_map = extract_band_to_2d(self.input_raster, left_band['band'])
        abs_map = extract_band_to_2d(self.input_raster, abs_band['band'])
        right_map = extract_band_to_2d(self.input_raster, right_band['band'])
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
