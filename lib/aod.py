#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AOD (Aerosol Optical Depth) estimation module for hyperspectral imagery.

This module provides functions to estimate Aerosol Optical Depth at 550nm
from hyperspectral data using various methods including the Dense Dark
Vegetation (DDV) method and other empirical approaches.

FIXES:
- Simplified band information extraction using r3.info metadata
- Optimized band extraction with proper 3D region handling
- Improved error handling and verbose output
"""

import os
import grass.script as gs
from wvc import get_band_info, find_nearest_band

# Wavelengths for common AOD estimation methods (in nm)
AOD_REFERENCE_WAVELENGTH = 550.0  # Standard reference wavelength for AOD
BLUE_BAND_WAVELENGTH = 490.0      # Typical blue band for AOD estimation
RED_BAND_WAVELENGTH = 670.0       # Typical red band for NDVI calculation
NIR_BAND_WAVELENGTH = 870.0       # Typical NIR band for NDVI calculation

# Default parameters for AOD estimation
DEFAULT_AOD_MIN = 0.01
DEFAULT_AOD_MAX = 1.5
DEFAULT_NDVI_THRESHOLD = 0.7      # NDVI threshold for DDV pixels
DEFAULT_DARK_PIXEL_THRESHOLD = 0.1 # Reflectance threshold for dark pixels


class AODEstimator:
    """Class for estimating Aerosol Optical Depth from hyperspectral data."""
    
    def __init__(self, input_raster, dem, sensor_config=None, verbose=False):
        """Initialize the AOD estimator.
        
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
        self.ndvi = None
        self.aod_map = None
        self.temp_maps = []  # Track temporary maps for cleanup
        
        # Get band information
        self._get_band_info()
    
    def _get_band_info(self):
        """Extract band information from the input raster metadata."""
        self.band_info = get_band_info(self.input_raster, verbose=self.verbose)
        
    def _find_nearest_band(self, target_wavelength):
        """Find the band closest to the target wavelength.
        
        Args:
            target_wavelength (float): Target wavelength in nm
            
        Returns:
            dict: Band information for the closest band
        """
        return find_nearest_band(self.band_info, target_wavelength)

    def extract_band_to_2d(self, band_num, output_map=None):
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
        """Calculate NDVI (Normalized Difference Vegetation Index)."""
        # Find nearest bands for red and NIR
        red_band = self._find_nearest_band(RED_BAND_WAVELENGTH)
        nir_band = self._find_nearest_band(NIR_BAND_WAVELENGTH)
        
        if self.verbose:
            gs.message(f"Using band {red_band['band']} ({red_band['wavelength']:.1f} nm) for red")
            gs.message(f"Using band {nir_band['band']} ({nir_band['wavelength']:.1f} nm) for NIR")
        
        # Extract bands from 3D raster
        red_map = self.extract_band_to_2d(red_band['band'], output_map='tmp_red')
        nir_map = self.extract_band_to_2d(nir_band['band'], output_map='tmp_nir')

        # Calculate NDVI: (NIR - Red) / (NIR + Red)
        ndvi_map = f"tmp_aod_ndvi_{os.getpid()}"
        expr = f"{ndvi_map} = float({nir_map} - {red_map}) / float({nir_map} + {red_map} + 0.0001)"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        
        self.ndvi = ndvi_map
        self.temp_maps.append(ndvi_map)
        self.temp_maps.append('temp_red')
        self.temp_maps.append('temp_nir')
        return ndvi_map
    
    def _estimate_aod_ddv(self):
        """Estimate AOD using the Dense Dark Vegetation (DDV) method."""
        if not self.ndvi:
            self._calculate_ndvi()
        
        # Find blue band for AOD estimation
        blue_band = self._find_nearest_band(BLUE_BAND_WAVELENGTH)
        
        if self.verbose:
            gs.message(f"Using band {blue_band['band']} ({blue_band['wavelength']:.1f} nm) for AOD estimation")
        
        # Extract blue band from 3D raster
        blue_map = self.extract_band_to_2d(blue_band['band'], output_map='temp_blue')
        self.temp_maps.append('temp_blue')
        
        # Create mask for DDV pixels (high NDVI and low blue reflectance)
        ddv_mask = f"tmp_aod_ddv_mask_{os.getpid()}"
        expr = (f"{ddv_mask} = if({self.ndvi} > {DEFAULT_NDVI_THRESHOLD} && "
                f"{blue_map} < {DEFAULT_DARK_PIXEL_THRESHOLD}, 1, null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(ddv_mask)
        
        # Count DDV pixels
        try:
            stats = gs.parse_command('r.univar', map=ddv_mask, flags='g')
            ddv_count = int(float(stats['n']))
        except:
            ddv_count = 0
        
        if ddv_count == 0:
            gs.warning("No suitable DDV pixels found for AOD estimation. Using default value.")
            return self._create_constant_aod(0.15)  # Default AOD
        
        if self.verbose:
            gs.message(f"Found {ddv_count} DDV pixels for AOD estimation")
        
        # For DDV pixels, estimate AOD from blue band reflectance
        aod_est = f"tmp_aod_est_{os.getpid()}"
        # Simplified DDV relationship: AOD ≈ -k * ln(rho_blue / rho_min)
        expr = f"{aod_est} = if({ddv_mask}, -0.15 * log(({blue_map} + 0.001) / 0.01), null())"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(aod_est)
        
        # Calculate mean AOD over DDV pixels
        try:
            stats = gs.parse_command('r.univar', map=aod_est, flags='g')
            mean_aod = float(stats['mean'])
        except:
            mean_aod = 0.15
        
        # Constrain to reasonable range
        mean_aod = max(DEFAULT_AOD_MIN, min(mean_aod, DEFAULT_AOD_MAX))
        
        if self.verbose:
            gs.message(f"Estimated AOD @ 550nm: {mean_aod:.3f}")
        
        # Create constant AOD map with the estimated value
        return self._create_constant_aod(mean_aod)
    
    def _create_constant_aod(self, value):
        """Create a constant AOD map with the given value."""
        aod_map = f"tmp_aod_550_{os.getpid()}"
        expr = f"{aod_map} = {value}"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(aod_map)
        return aod_map
    
    def estimate_aod(self, method='ddv'):
        """Estimate AOD using the specified method.
        
        Args:
            method (str): AOD estimation method ('ddv', 'dark_pixel', 'constant')
            
        Returns:
            tuple: (aod_map, mean_aod) where aod_map is the raster name and mean_aod is the value
        """
        if method == 'ddv':
            self.aod_map = self._estimate_aod_ddv()
        elif method == 'dark_pixel':
            # Similar to DDV but with different thresholds and no vegetation requirement
            gs.warning("Dark pixel method not fully implemented, using DDV method")
            self.aod_map = self._estimate_aod_ddv()
        elif method == 'constant':
            default_aod = 0.15
            self.aod_map = self._create_constant_aod(default_aod)
        else:
            raise ValueError(f"Unknown AOD estimation method: {method}")
        
        # Get mean AOD value
        try:
            stats = gs.parse_command('r.univar', map=self.aod_map, flags='g')
            mean_aod = float(stats['mean'])
        except:
            mean_aod = 0.15
        
        if self.verbose:
            gs.message(f"AOD estimation complete. Mean AOD: {mean_aod:.3f}")
        
        return self.aod_map, mean_aod
    
    def cleanup(self):
        """Clean up temporary maps."""
        if self.verbose:
            gs.message("Cleaning up temporary maps...")
        
        for map_name in self.temp_maps:
            if gs.find_file(map_name, element='cell')['file']:
                gs.run_command('g.remove', flags='f', type='raster', 
                             name=map_name, quiet=True)


def estimate_aod(input_raster, dem, method='ddv', sensor_config=None, verbose=False):
    """Convenience function to estimate AOD from hyperspectral data.
    
    Args:
        input_raster (str): Name of the input 3D hyperspectral raster
        dem (str): Name of the digital elevation model (DEM) raster
        method (str, optional): AOD estimation method. Defaults to 'ddv'.
        sensor_config (dict, optional): Sensor configuration dictionary
        verbose (bool, optional): Enable verbose output
        
    Returns:
        tuple: (aod_map, mean_aod) where aod_map is the raster name and mean_aod is the value
    """
    estimator = AODEstimator(input_raster, dem, sensor_config, verbose)
    try:
        aod_map, mean_aod = estimator.estimate_aod(method)
        
        # Create a permanent copy
        permanent_aod = f"{input_raster}_aod_550"
        gs.run_command('g.copy', raster=f"{aod_map},{permanent_aod}", 
                      overwrite=True, quiet=not verbose)
        
        return permanent_aod, mean_aod
    finally:
        estimator.cleanup()


if __name__ == "__main__":
    # Example usage
    import argparse
    
    parser = argparse.ArgumentParser(description='Estimate AOD from hyperspectral data')
    parser.add_argument('input', help='Input 3D hyperspectral raster')
    parser.add_argument('dem', help='Digital Elevation Model (DEM) raster')
    parser.add_argument('--method', default='ddv', 
                       choices=['ddv', 'dark_pixel', 'constant'],
                       help='AOD estimation method')
    parser.add_argument('--output', help='Output AOD map name')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    aod_map, mean_aod = estimate_aod(args.input, args.dem, args.method, None, args.verbose)
    
    if args.output and args.output != aod_map:
        gs.run_command('g.rename', raster=f"{aod_map},{args.output}", overwrite=True)
        aod_map = args.output
    
    print(f"AOD map created: {aod_map}")
    print(f"Mean AOD @ 550nm: {mean_aod:.3f}")
