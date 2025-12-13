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
            
            # Set the 3D region back
            gs.run_command('g.region', raster_3d=self.input_raster, quiet=True)

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
        """Estimate AOD using the Dense Dark Vegetation (DDV) method with fallback to spectral shape."""
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
        
        if ddv_count < 100:  # Minimum number of DDV pixels needed
            if self.verbose:
                gs.message(f"Only {ddv_count} DDV pixels found. Falling back to spectral shape method.")
            return self._estimate_aod_spectral_shape()
        
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
            if self.verbose:
                gs.message(f"DDV AOD estimation complete. Mean AOD: {mean_aod:.3f}")
            return aod_est, mean_aod
        except:
            gs.warning("Error calculating AOD statistics. Falling back to spectral shape method.")
            return self._estimate_aod_spectral_shape()  
    
    def _estimate_aod_spectral_shape(self):
        """Estimate AOD using spectral shape analysis when DDV pixels are not available.
        
        This method uses the spectral shape in the NIR-SWIR region where the surface 
        reflectance is less affected by aerosols, and the atmospheric path radiance 
        is minimal.
        """
        if self.verbose:
            gs.message("Using spectral shape method for AOD estimation")
        
        try:
            # Find key wavelengths for AOD estimation
            blue_band = self._find_nearest_band(470)  # Blue band for AOD sensitivity
            green_band = self._find_nearest_band(550)  # Green band
            red_band = self._find_nearest_band(670)   # Red band
            nir_band = self._find_nearest_band(870)   # NIR band
            swir1_band = self._find_nearest_band(1640)  # SWIR1 band
            swir2_band = self._find_nearest_band(2130)  # SWIR2 band
            
            # Extract required bands
            blue_map = self.extract_band_to_2d(blue_band['band'], output_map='temp_blue')
            green_map = self.extract_band_to_2d(green_band['band'], output_map='temp_green')
            red_map = self.extract_band_to_2d(red_band['band'], output_map='temp_red')
            nir_map = self.extract_band_to_2d(nir_band['band'], output_map='temp_nir')
            swir1_map = self.extract_band_to_2d(swir1_band['band'], output_map='temp_swir1') if swir1_band else None
            swir2_map = self.extract_band_to_2d(swir2_band['band'], output_map='temp_swir2') if swir2_band else None
            
            # Calculate spectral indices
            # 1. Aerosol Free Vegetation Index (AFRI) - less sensitive to aerosols
            afri_map = f"tmp_afri_{os.getpid()}"
            expr = f"{afri_map} = float({nir_map} - 0.5 * {swir1_map}) / float({nir_map} + 0.5 * {swir1_map} + 0.0001)" if swir1_map else None
            if expr:
                gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
                self.temp_maps.append(afri_map)
            
            # 2. Normalized Difference Vegetation Index (NDVI)
            ndvi_map = f"tmp_ndvi_{os.getpid()}"
            expr = f"{ndvi_map} = float({nir_map} - {red_map}) / float({nir_map} + {red_map} + 0.0001)"
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(ndvi_map)
            
            # 3. Blue/Red ratio - sensitive to aerosol content
            blue_red_ratio = f"tmp_brr_{os.getpid()}"
            expr = f"{blue_red_ratio} = float({blue_map} + 0.01) / float({red_map} + 0.01)"
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(blue_red_ratio)
            
            # Create a mask for valid pixels (non-water, non-cloud)
            valid_mask = f"tmp_valid_{os.getpid()}"
            if swir1_map and swir2_map:
                # Use SWIR-based cloud/water mask
                expr = (f"{valid_mask} = if({swir1_map} > 0.05 && {swir2_map} > 0.05 && "
                       f"{ndvi_map} > -0.1 && {ndvi_map} < 0.8, 1, null())")
            else:
                # Fallback to NDVI-based mask
                expr = f"{valid_mask} = if({ndvi_map} > -0.2 && {ndvi_map} < 0.8, 1, null())"
            
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(valid_mask)
            
            # Estimate AOD using an empirical relationship based on blue/red ratio
            aod_est = f"tmp_aod_est_{os.getpid()}"
            expr = (f"{aod_est} = if({valid_mask}, "
                   f"0.1 + 0.8 * (1.0 - exp(-2.5 * {blue_red_ratio})), null())")
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(aod_est)
            
            # Calculate mean AOD over valid pixels
            try:
                stats = gs.parse_command('r.univar', map=aod_est, flags='g')
                mean_aod = float(stats['mean'])
                if self.verbose:
                    gs.message(f"Spectral shape AOD estimation complete. Mean AOD: {mean_aod:.3f}")
                return aod_est, mean_aod
            except:
                gs.warning("Error calculating AOD statistics. Using default value.")
                return self._create_constant_aod(0.15)
                
        except Exception as e:
            gs.warning(f"Error in spectral shape AOD estimation: {str(e)}. Using default value.")
            return self._create_constant_aod(0.15)
    
    def _create_constant_aod(self, value):
        """Create a constant AOD map with the given value."""
        aod_map = f"tmp_aod_550_{os.getpid()}"
        expr = f"{aod_map} = {value}"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(aod_map)
        return aod_map
    
    def estimate_aod(self, method='auto'):
        """Estimate AOD using the specified method.
        
        Args:
            method (str): AOD estimation method:
                - 'auto': Automatically select the best available method
                - 'ddv': Dense Dark Vegetation method (default)
                - 'spectral': Spectral shape analysis
                - 'dark_pixel': Dark pixel method
                - 'constant': Use a constant AOD value
                
        Returns:
            tuple: (aod_map, mean_aod) where aod_map is the raster name and mean_aod is the value
        """
        if method == 'auto':
            # First try DDV method
            try:
                aod_map, mean_aod = self._estimate_aod_ddv()
                if mean_aod > 0:  # If we got a valid result
                    return aod_map, mean_aod
            except Exception as e:
                if self.verbose:
                    gs.message(f"DDV method failed: {str(e)}. Trying spectral method...")
            
            # Fall back to spectral method
            try:
                return self._estimate_aod_spectral_shape()
            except Exception as e:
                if self.verbose:
                    gs.message(f"Spectral method failed: {str(e)}. Using constant AOD...")
                return self._create_constant_aod(0.15), 0.15
                
        elif method == 'ddv':
            return self._estimate_aod_ddv()
            
        elif method == 'spectral':
            return self._estimate_aod_spectral_shape()
            
        elif method == 'dark_pixel':
            return self._estimate_aod_dark_pixel()
            
        elif method == 'constant':
            default_aod = 0.15
            aod_map = self._create_constant_aod(default_aod)
            return aod_map, default_aod
            
        else:
            raise ValueError(f"Unknown AOD estimation method: {method}")
    
    def _estimate_aod_dark_pixel(self):
        """Estimate AOD using dark pixel method (similar to DDV but no vegetation requirement)."""
        if self.verbose:
            gs.message("Using dark pixel method for AOD estimation")
        
        try:
            # Find blue band for AOD estimation
            blue_band = self._find_nearest_band(470)  # Blue band for AOD sensitivity
            swir2_band = self._find_nearest_band(2200)  # SWIR2 band for dark pixel detection
            
            # Extract bands
            blue_map = self.extract_band_to_2d(blue_band['band'], output_map='temp_blue')
            swir2_map = self.extract_band_to_2d(swir2_band['band'], output_map='temp_swir2') if swir2_band else None
            
            # Create mask for dark pixels (low SWIR2 reflectance)
            dark_mask = f"tmp_dark_mask_{os.getpid()}"
            if swir2_map:
                expr = f"{dark_mask} = if({swir2_map} < 0.1, 1, null())"  # Very dark pixels in SWIR
            else:
                # Fallback to using blue band if SWIR is not available
                expr = f"{dark_mask} = if({blue_map} < 0.05, 1, null())"  # Very dark pixels in blue
            
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(dark_mask)
            
            # Count dark pixels
            try:
                stats = gs.parse_command('r.univar', map=dark_mask, flags='g')
                dark_count = int(float(stats['n']))
            except:
                dark_count = 0
            
            if dark_count < 100:  # Minimum number of dark pixels needed
                if self.verbose:
                    gs.message(f"Only {dark_count} dark pixels found. Falling back to spectral method.")
                return self._estimate_aod_spectral_shape()
            
            if self.verbose:
                gs.message(f"Found {dark_count} dark pixels for AOD estimation")
            
            # For dark pixels, estimate AOD from blue band reflectance
            aod_est = f"tmp_aod_est_{os.getpid()}"
            # Empirical relationship for dark pixels
            expr = f"{aod_est} = if({dark_mask}, 0.1 + 0.5 * (1.0 - exp(-5.0 * {blue_map})), null())"
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=not self.verbose)
            self.temp_maps.append(aod_est)
            
            # Calculate mean AOD over dark pixels
            stats = gs.parse_command('r.univar', map=aod_est, flags='g')
            mean_aod = float(stats['mean'])
            if self.verbose:
                gs.message(f"Dark pixel AOD estimation complete. Mean AOD: {mean_aod:.3f}")
            return aod_est, mean_aod
            
        except Exception as e:
            if self.verbose:
                gs.warning(f"Error in dark pixel AOD estimation: {str(e)}. Falling back to constant AOD.")
            return self._create_constant_aod(0.15), 0.15
    
    def cleanup(self):
        """Clean up temporary maps."""
        if self.verbose:
            gs.message("Cleaning up temporary maps...")
        
        for map_name in self.temp_maps:
            if gs.find_file(map_name, element='cell')['file']:
                gs.run_command('g.remove', flags='f', type='raster', 
                             name=map_name, quiet=True)


def estimate_aod(input_raster, dem, method='auto', sensor_config=None, verbose=False):
    """Convenience function to estimate AOD from hyperspectral data.
    
    Args:
        input_raster (str): Name of the input 3D hyperspectral raster
        dem (str): Name of the digital elevation model (DEM) raster
        method (str, optional): AOD estimation method. Defaults to 'auto'.
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
    parser.add_argument('--method', default='auto', 
                       choices=['ddv', 'dark_pixel', 'constant', 'spectral'],
                       help='AOD estimation method')
    parser.add_argument('--output', help='Output AOD map name')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    aod_map, mean_aod = estimate_aod(args.input, args.dem, args.method, None, args.verbose)
    
    if args.output and args.output != aod_map:
        gs.run_command('g.rename', raster=f"{aod_map},{args.output}", overwrite=True)
        aod_map = args.output
    
    if gs.verbosity() > 0:
        gs.message(f"AOD map created: {aod_map}")
        gs.message(f"Mean AOD @ 550nm: {mean_aod:.3f}")
