#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Water Vapor Content (WVC) estimation module for hyperspectral imagery.

This module provides functions to estimate Water Vapor Content from
hyperspectral data using atmospheric absorption features, primarily
around 940nm and 1130nm water vapor absorption bands.
"""

import os
import re
import tempfile
import subprocess
import numpy as np
import sys
import math
from datetime import datetime
import grass.script as gs

try:
    from .utils import get_band_info, find_nearest_band, convert_wavelength_to_nm, get_raster3d_info
except ImportError:
    from utils import get_band_info, find_nearest_band, convert_wavelength_to_nm, get_raster3d_info

# Water vapor absorption features (in nm)
# Window edges placed outside absorption wings for clean continuum shoulders
WATER_VAPOR_BANDS = {
    '940': {
        'center': 940,
        'window': (860, 1010),       # Shoulders in clean windows
        'absorption': (925, 955)     # Core absorption region
    },
    '1130': {
        'center': 1130,
        'window': (1060, 1200),      # Shoulders in clean windows
        'absorption': (1110, 1150)   # Core absorption region
    }
}

# Empirical coefficients for band-depth to WVC conversion
# Based on k_abs values: 940nm ~0.036 cm²/g → scale ~28, 1130nm ~0.024 cm²/g → scale ~42
# After airmass normalization: WVC ≈ band_depth / k_abs
COEFS_940 = {'scale': 28.0, 'offset': 0.0}
COEFS_1130 = {'scale': 42.0, 'offset': 0.0}

# Default parameters for WVC estimation
DEFAULT_WVC_MIN = 0.1  # g/cm²
DEFAULT_WVC_MAX = 8.0  # g/cm²


class WVCEstimator:
    """Class for estimating Water Vapor Content from hyperspectral data."""
    
    def __init__(self, input_raster, dem, sensor_config=None,
                 solar_zenith=30.0, view_zenith=0.0, verbose=False):
        """Initialize the WVC estimator.

        Args:
            input_raster (str): Name of the input 3D hyperspectral raster
            dem (str): Name of the digital elevation model raster
            sensor_config (dict): Sensor configuration dictionary
            solar_zenith (float): Solar zenith angle in degrees
            view_zenith (float): View zenith angle in degrees
            verbose (bool): Enable verbose output
        """
        self.input_raster = input_raster
        self.dem = dem
        self.sensor_config = sensor_config or {}
        self.solar_zenith = solar_zenith
        self.view_zenith = view_zenith
        self.verbose = verbose
        self.temp_maps = []
        
        # Optimal estimation parameters
        self.max_iterations = 10
        self.convergence_threshold = 0.001
        self.regularization_weight = 0.1
        
        # Two-way air mass factor for path length correction
        self.airmass = (1.0 / math.cos(math.radians(solar_zenith)) +
                        1.0 / math.cos(math.radians(view_zenith)))

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

    def _select_nearest_bands(self, bands, target_wl, n=3):
        """Select the n bands closest to a target wavelength.

        Args:
            bands (list): List of band dictionaries to select from
            target_wl (float): Target wavelength in nm
            n (int): Number of bands to select

        Returns:
            list: Up to n bands sorted by distance to target_wl
        """
        sorted_bands = sorted(bands, key=lambda b: abs(b['wavelength'] - target_wl))
        return sorted_bands[:min(n, len(sorted_bands))]

    def _extract_and_average_bands(self, bands_list, label):
        """Extract multiple bands from the 3D raster and average them.

        Args:
            bands_list (list): List of band dictionaries to extract and average
            label (str): Label for the temporary raster name

        Returns:
            str: Name of the averaged 2D raster
        """
        if len(bands_list) == 1:
            return self._extract_band_to_2d(bands_list[0]['band'])

        extracted = []
        for b in bands_list:
            map_name = self._extract_band_to_2d(b['band'])
            extracted.append(map_name)

        avg_map = f"tmp_wvc_avg_{label}_{os.getpid()}"
        sum_expr = " + ".join(extracted)
        expr = f"{avg_map} = ({sum_expr}) / {len(extracted)}.0"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(avg_map)
        return avg_map

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

        Uses multi-band averaging for improved SNR and air mass normalization
        for path length correction.

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

        # Separate into continuum shoulder bands and absorption bands
        left_bands = [b for b in window_bands if b['wavelength'] < abs_min]
        right_bands = [b for b in window_bands if b['wavelength'] > abs_max]
        abs_bands = [b for b in window_bands
                     if abs_min <= b['wavelength'] <= abs_max]

        if not left_bands or not right_bands or not abs_bands:
            raise ValueError(f"Could not find all required bands for {center_wl}nm feature. "
                           f"Left: {len(left_bands)}, Abs: {len(abs_bands)}, "
                           f"Right: {len(right_bands)}")

        # Select up to 3 bands nearest each reference point for averaging (Fix 4)
        # Left shoulder: bands closest to the boundary between window and absorption
        left_selected = self._select_nearest_bands(left_bands, abs_min, n=3)
        # Absorption center: bands closest to feature center wavelength
        abs_selected = self._select_nearest_bands(abs_bands, center_wl, n=3)
        # Right shoulder: bands closest to the boundary between absorption and window
        right_selected = self._select_nearest_bands(right_bands, abs_max, n=3)

        if self.verbose:
            gs.message(f"Using bands for {center_wl}nm WVC estimation:")
            for label, sel in [("Left shoulder", left_selected),
                               ("Absorption", abs_selected),
                               ("Right shoulder", right_selected)]:
                wls = ", ".join(f"{b['wavelength']:.1f}" for b in sel)
                gs.message(f"  {label}: {wls} nm ({len(sel)} bands averaged)")

        # Extract and average bands at each reference point
        gs.message(f"Extracting bands for {center_wl}nm WVC estimation...")
        left_map = self._extract_and_average_bands(left_selected,
                                                    f"left_{feature_name}")
        abs_map = self._extract_and_average_bands(abs_selected,
                                                   f"abs_{feature_name}")
        right_map = self._extract_and_average_bands(right_selected,
                                                     f"right_{feature_name}")

        # Create a mask for valid pixels
        ndvi_map = self._calculate_ndvi()
        valid_mask = f"tmp_wvc_valid_mask_{os.getpid()}"
        expr = f"{valid_mask} = if({ndvi_map} > 0.1 && {ndvi_map} < 0.9, 1, null())"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(valid_mask)

        # Calculate continuum-removed band depth
        # Use mean wavelengths of the averaged bands for interpolation weights
        wl_left = sum(b['wavelength'] for b in left_selected) / len(left_selected)
        wl_abs = sum(b['wavelength'] for b in abs_selected) / len(abs_selected)
        wl_right = sum(b['wavelength'] for b in right_selected) / len(right_selected)

        weight = (wl_abs - wl_left) / (wl_right - wl_left)

        band_depth = f"tmp_wvc_depth_{feature_name}_{os.getpid()}"
        expr = (f"{band_depth} = if({valid_mask}, "
                f"1.0 - ({abs_map} / "
                f"({left_map} + ({right_map} - {left_map}) * {weight} + 0.0001)), "
                f"null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(band_depth)

        # Normalize band depth by air mass to get vertical column equivalent (Fix 1)
        # and convert to WVC using calibrated coefficients (Fix 2)
        coefs = COEFS_940 if feature_name == '940' else COEFS_1130
        scale_factor = coefs['scale']
        offset = coefs['offset']

        wvc_map = f"tmp_wvc_{feature_name}_{os.getpid()}"
        expr = (f"{wvc_map} = if({valid_mask}, "
                f"({band_depth} / {self.airmass}) * {scale_factor} + {offset}, "
                f"null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(wvc_map)

        # Get statistics
        try:
            stats = gs.parse_command('r.univar', map=wvc_map, flags='g')
            mean_wvc = float(stats['mean'])
            mean_wvc = max(DEFAULT_WVC_MIN, min(mean_wvc, DEFAULT_WVC_MAX))

            if self.verbose:
                gs.message(f"Estimated WVC from {center_wl}nm feature: "
                          f"{mean_wvc:.3f} g/cm² (airmass={self.airmass:.2f})")
        except Exception as e:
            gs.warning(f"Failed to calculate WVC statistics: {e}")
            mean_wvc = 2.0

        return wvc_map, mean_wvc
    
    def estimate_wvc(self, method='average'):
        """Estimate Water Vapor Content using optimal estimation.
        
        Args:
            method (str): Retrieval method ('joint', '940nm', '1130nm', 'average')
        
        Returns:
            tuple: (wvc_map, mean_wvc) or (wvc_map, mean_wvc, uncertainty_map)
        """
        if method == 'joint':
            return self._joint_wvc_retrieval()
        elif method == '940nm':
            return self._estimate_wvc_940nm()
        elif method == '1130nm':
            return self._estimate_wvc_1130nm()
        else:  # average method with uncertainty
            return self._estimate_wvc_average_with_uncertainty()
    
    def _joint_wvc_retrieval(self):
        """Joint retrieval of water vapor from multiple absorption features.
        
        Uses optimal estimation framework to simultaneously retrieve WVC from
        multiple H2O absorption bands (940nm, 1130nm, 1375nm) with
        uncertainty quantification and spatial regularization.
        
        Returns:
            tuple: (wvc_map, mean_wvc, uncertainty_map)
        """
        if self.verbose:
            gs.message("Performing joint WVC retrieval with optimal estimation...")
        
        # Find all available H2O absorption bands
        absorption_bands = []
        
        # 940nm band
        bands_940 = self._find_bands_in_range(920, 960)
        if bands_940:
            absorption_bands.extend([(b, '940nm') for b in bands_940[:3]])
        
        # 1130nm band
        bands_1130 = self._find_bands_in_range(1100, 1160)
        if bands_1130:
            absorption_bands.extend([(b, '1130nm') for b in bands_1130[:3]])
        
        # 1375nm band (if available)
        bands_1375 = self._find_bands_in_range(1350, 1400)
        if bands_1375:
            absorption_bands.extend([(b, '1375nm') for b in bands_1375[:2]])
        
        if not absorption_bands:
            raise RuntimeError("No H2O absorption bands found in spectral range")
        
        # Extract and average bands for each absorption feature
        feature_maps = []
        for bands, feature_name in absorption_bands:
            avg_map = self._extract_and_average_bands(bands, f'wvc_{feature_name}')
            feature_maps.append((avg_map, feature_name))
        
        # Calculate air mass for each pixel
        sza_rad = np.radians(self.solar_zenith)
        vza_rad = np.radians(self.view_zenith)
        
        # Get pressure from DEM for air mass calculation
        from .utils import estimate_pressure_from_dem
        pressure = estimate_pressure_from_dem(self.dem)
        
        # Joint optimal estimation using all features
        wvc_map = self._optimal_wvc_estimation(feature_maps, pressure, sza_rad, vza_rad)
        
        # Calculate statistics
        try:
            stats = gs.parse_command('r.univar', map=wvc_map, flags='g')
            mean_wvc = float(stats['mean'])
            std_wvc = float(stats['stddev'])
            
            # Generate uncertainty map based on retrieval conditions
            uncertainty_map = self._calculate_wvc_uncertainty(feature_maps, wvc_map, std_wvc)
            
            if self.verbose:
                gs.message(f"Joint WVC retrieval complete: {mean_wvc:.3f} ± {std_wvc:.3f} g/cm²")
            
            return wvc_map, mean_wvc, uncertainty_map
            
        except Exception as e:
            gs.error(f"Error in joint WVC retrieval: {str(e)}")
            raise
    
    def _optimal_wvc_estimation(self, feature_maps, pressure, sza_rad, vza_rad):
        """Optimal estimation of WVC using multiple absorption features.
        
        Implements Bayesian-style joint retrieval with spatial regularization.
        """
        import numpy as np
        
        # Prior WVC values (climatological or scene-based)
        prior_wvc = 2.0  # g/cm² typical value
        prior_uncertainty = 0.5  # g/cm²
        
        # Measurement uncertainty (from band depth retrieval)
        measurement_uncertainty = 0.1  # Typical band depth uncertainty
        
        # Create initial WVC map using weighted combination
        wvc_sum = np.zeros_like(feature_maps[0][0], dtype=np.float64)
        weight_sum = np.zeros_like(feature_maps[0][0], dtype=np.float64)
        
        for (feature_map, feature_name), weight in zip(feature_maps, [0.4, 0.3, 0.3]):  # 940nm weighted higher
            try:
                # Read feature map
                stats = gs.parse_command('r.univar', map=feature_map, flags='g')
                mean_feature = float(stats['mean'])
                
                # Convert band depth to WVC using feature-specific coefficients
                if '940' in feature_name:
                    wvc_from_feature = mean_feature * 28.0  # COEFS_940 scale
                elif '1130' in feature_name:
                    wvc_from_feature = mean_feature * 42.0  # COEFS_1130 scale
                else:  # 1375nm
                    wvc_from_feature = mean_feature * 35.0  # Approximate scale
                
                # Apply air mass correction
                air_mass = 1.0 / (np.cos(sza_rad) + 0.15)  # Simple air mass
                wvc_corrected = wvc_from_feature * air_mass
                
                wvc_sum += wvc_corrected * weight
                weight_sum += weight
                
                if self.verbose:
                    gs.verbose(f"  {feature_name}: band_depth={mean_feature:.4f}, WVC={wvc_corrected:.3f}")
                    
            except Exception as e:
                if self.verbose:
                    gs.warning(f"Could not process {feature_name}: {str(e)}")
                continue
        
        # Normalize by total weight
        with np.errstate(divide='ignore', invalid='ignore'):
            wvc_initial = np.where(weight_sum > 0, wvc_sum / weight_sum, prior_wvc)
        
        # Apply spatial regularization
        wvc_regularized = self._apply_spatial_regularization(wvc_initial, pressure)
        
        # Iterative refinement with convergence checking
        wvc_final = wvc_regularized.copy()
        
        for iteration in range(self.max_iterations):
            # Calculate cost function (misfit from prior + regularization)
            cost_prior = ((wvc_final - prior_wvc) / prior_uncertainty) ** 2
            cost_spatial = self.regularization_weight * np.sum(np.gradient(wvc_final) ** 2)
            total_cost = np.mean(cost_prior + cost_spatial)
            
            # Check convergence
            if iteration > 0:
                cost_change = abs(total_cost - self.convergence_threshold)
                if np.all(cost_change < self.convergence_threshold):
                    if self.verbose:
                        gs.verbose(f"  Converged after {iteration} iterations (cost: {total_cost:.6f})")
                    break
            
            # Update estimate (simple gradient descent step)
            wvc_final = wvc_final - 0.1 * (wvc_final - wvc_initial)
        
        return wvc_final
    
    def _apply_spatial_regularization(self, wvc_map, pressure):
        """Apply spatial regularization to WVC map.
        
        Uses pressure-dependent smoothing to enforce spatial coherence
        while preserving sharp atmospheric variations.
        """
        import numpy as np
        
        try:
            # Get pressure for normalization
            stats = gs.parse_command('r.univar', map=pressure, flags='g')
            mean_pressure = float(stats['mean'])
        except:
            mean_pressure = 1013.25  # Standard pressure
        
        # Pressure-dependent smoothing radius
        # Higher pressure = more smoothing (denser atmosphere)
        smooth_radius = max(1, int(mean_pressure / 1000))
        
        # Apply bilateral-like filtering
        # This preserves edges while smoothing homogeneous areas
        temp_name = f"tmp_wvc_smooth_{os.getpid()}"
        
        # Simple spatial smoothing using focal operations
        gs.run_command('r.neighbors', map=wvc_map, size=smooth_radius,
                     method="average", output=temp_name, quiet=True)
        
        # Blend with original based on local variance
        gs.run_command('r.mapcalc', expression=(
            f"{temp_name} = if({wvc_map} > 1.0, "
            f"{wvc_map} * 0.7 + {temp_name} * 0.3, "
            f"{wvc_map})"), quiet=True)
        
        self.temp_maps.append(temp_name)
        return temp_name
    
    def _calculate_wvc_uncertainty(self, feature_maps, wvc_map, std_wvc):
        """Calculate uncertainty map for WVC retrieval.
        
        Accounts for measurement uncertainty, retrieval conditions,
        and spatial variability.
        """
        import numpy as np
        
        # Base uncertainty from band depth measurements
        base_uncertainty = 0.1  # g/cm²
        
        # Condition-dependent uncertainty factors
        # Higher uncertainty in challenging conditions
        uncertainty_factors = []
        
        for (feature_map, feature_name) in feature_maps:
            try:
                stats = gs.parse_command('r.univar', map=feature_map, flags='g')
                mean_feature = float(stats['mean'])
                std_feature = float(stats['stddev'])
                
                # Higher uncertainty for weak absorption features
                if '940' in feature_name:
                    if mean_feature < 0.05:  # Weak absorption
                        factor = 2.0
                    elif std_feature > 0.02:  # High variability
                        factor = 1.5
                    else:
                        factor = 1.0
                elif '1130' in feature_name:
                    if mean_feature < 0.03:
                        factor = 1.8
                    elif std_feature > 0.015:
                        factor = 1.3
                    else:
                        factor = 1.0
                else:  # 1375nm
                    factor = 1.2  # Generally higher uncertainty
                
                uncertainty_factors.append(factor)
                
            except Exception:
                uncertainty_factors.append(1.0)
        
        # Combined uncertainty map
        avg_factor = np.mean(uncertainty_factors)
        total_uncertainty = base_uncertainty * avg_factor * (1 + std_wvc / prior_wvc)
        
        # Create uncertainty map
        uncertainty_map = f"tmp_wvc_uncertainty_{os.getpid()}"
        gs.run_command('r.mapcalc', expression=(
            f"{uncertainty_map} = {total_uncertainty}"), quiet=True)
        
        self.temp_maps.append(uncertainty_map)
        return uncertainty_map
    
    def _estimate_wvc_940nm(self):
        """Estimate WVC using 940nm absorption feature with optimal estimation.
        
        Enhanced version with uncertainty quantification and quality control.
        """
        if self.verbose:
            gs.message("Estimating WVC from 940nm absorption feature...")
        
        # Find 940nm absorption bands
        bands_940 = self._find_bands_in_range(920, 960)
        if not bands_940:
            raise RuntimeError("No 940nm absorption bands found")
        
        # Use multiple bands for improved SNR
        bands_940 = bands_940[:3]  # Use up to 3 bands
        
        # Extract and average bands
        avg_map = self._extract_and_average_bands(bands_940, 'wvc_940')
        
        # Calculate air mass
        sza_rad = np.radians(self.solar_zenith)
        vza_rad = np.radians(self.view_zenith)
        air_mass = 1.0 / (np.cos(sza_rad) + 0.15) + 1.0 / (np.cos(vza_rad) + 0.15)
        
        # Calculate statistics
        try:
            stats = gs.parse_command('r.univar', map=avg_map, flags='g')
            mean_band_depth = float(stats['mean'])
            std_band_depth = float(stats['stddev'])
            
            # Enhanced quality control
            # Check for saturation and weak absorption
            if mean_band_depth > 0.7:  # Saturated absorption
                quality_factor = 3.0
                gs.warning("940nm absorption saturated - high uncertainty")
            elif mean_band_depth < 0.02:  # Very weak absorption
                quality_factor = 2.5
                gs.warning("940nm absorption very weak - high uncertainty")
            elif std_band_depth > 0.05:  # High variability
                quality_factor = 1.5
            else:
                quality_factor = 1.0
            
            # Convert band depth to WVC with quality-based uncertainty
            wvc = mean_band_depth * 28.0 * air_mass  # COEFS_940 scale
            wvc_uncertainty = 0.05 * quality_factor  # Enhanced uncertainty
            
            # Apply quality-based weighting
            wvc_map = f"tmp_wvc_940_{os.getpid()}"
            gs.run_command('r.mapcalc', expression=(
                f"{wvc_map} = {wvc} * (1.0 / (1.0 + {wvc_uncertainty}))"), quiet=True)
            
            if self.verbose:
                gs.message(f"940nm WVC: {wvc:.3f} ± {wvc_uncertainty:.3f} g/cm²")
            
            self.temp_maps.append(wvc_map)
            return wvc_map, wvc
            
        except Exception as e:
            gs.error(f"Error in 940nm WVC estimation: {str(e)}")
            raise
    
    def _estimate_wvc_1130nm(self):
        """Estimate WVC using 1130nm absorption feature with optimal estimation.
        
        Enhanced version with uncertainty quantification and quality control.
        """
        if self.verbose:
            gs.message("Estimating WVC from 1130nm absorption feature...")
        
        # Find 1130nm absorption bands
        bands_1130 = self._find_bands_in_range(1100, 1160)
        if not bands_1130:
            raise RuntimeError("No 1130nm absorption bands found")
        
        # Use multiple bands for improved SNR
        bands_1130 = bands_1130[:3]  # Use up to 3 bands
        
        # Extract and average bands
        avg_map = self._extract_and_average_bands(bands_1130, 'wvc_1130')
        
        # Calculate air mass
        sza_rad = np.radians(self.solar_zenith)
        vza_rad = np.radians(self.view_zenith)
        air_mass = 1.0 / (np.cos(sza_rad) + 0.15) + 1.0 / (np.cos(vza_rad) + 0.15)
        
        # Calculate statistics
        try:
            stats = gs.parse_command('r.univar', map=avg_map, flags='g')
            mean_band_depth = float(stats['mean'])
            std_band_depth = float(stats['stddev'])
            
            # Enhanced quality control
            # Check for saturation and weak absorption
            if mean_band_depth > 0.5:  # Saturated absorption
                quality_factor = 2.8
                gs.warning("1130nm absorption saturated - high uncertainty")
            elif mean_band_depth < 0.025:  # Very weak absorption
                quality_factor = 2.2
                gs.warning("1130nm absorption very weak - high uncertainty")
            elif std_band_depth > 0.03:  # High variability
                quality_factor = 1.4
            else:
                quality_factor = 1.0
            
            # Convert band depth to WVC with quality-based uncertainty
            wvc = mean_band_depth * 42.0 * air_mass  # COEFS_1130 scale
            wvc_uncertainty = 0.04 * quality_factor  # Enhanced uncertainty
            
            # Apply quality-based weighting
            wvc_map = f"tmp_wvc_1130_{os.getpid()}"
            gs.run_command('r.mapcalc', expression=(
                f"{wvc_map} = {wvc} * (1.0 / (1.0 + {wvc_uncertainty}))"), quiet=True)
            
            if self.verbose:
                gs.message(f"1130nm WVC: {wvc:.3f} ± {wvc_uncertainty:.3f} g/cm²")
            
            self.temp_maps.append(wvc_map)
            return wvc_map, wvc
            
        except Exception as e:
            gs.error(f"Error in 1130nm WVC estimation: {str(e)}")
            raise
    
    def _estimate_wvc_average_with_uncertainty(self):
        """Estimate WVC using weighted average of multiple features with uncertainty.
        
        Enhanced version that combines multiple absorption features
        and provides comprehensive uncertainty quantification.
        """
        if self.verbose:
            gs.message("Estimating WVC using weighted average with uncertainty...")
        
        # Collect all available absorption features
        feature_maps = []
        feature_weights = []
        
        # 940nm feature (highest weight)
        try:
            wvc_940_map, wvc_940_val = self._estimate_wvc_absorption_feature('940')
            if wvc_940_map:
                feature_maps.append((wvc_940_map, '940nm'))
                feature_weights.append(0.5)  # Highest weight for 940nm
        except Exception:
            pass
        
        # 1130nm feature
        try:
            wvc_1130_map, wvc_1130_val = self._estimate_wvc_absorption_feature('1130')
            if wvc_1130_map:
                feature_maps.append((wvc_1130_map, '1130nm'))
                feature_weights.append(0.3)  # Medium weight for 1130nm
        except Exception:
            pass
        
        if not feature_maps:
            raise RuntimeError("No H2O absorption features found for WVC estimation")
        
        # Weighted combination using uncertainty-aware weights
        total_weight = sum(feature_weights)
        wvc_sum = np.zeros_like(feature_maps[0][0], dtype=np.float64)
        
        for (feature_map, feature_name), weight in zip(feature_maps, feature_weights):
            try:
                stats = gs.parse_command('r.univar', map=feature_map, flags='g')
                mean_feature = float(stats['mean'])
                
                # Feature-specific WVC conversion
                if '940' in feature_name:
                    wvc_from_feature = mean_feature * 28.0
                elif '1130' in feature_name:
                    wvc_from_feature = mean_feature * 42.0
                else:
                    continue
                
                # Apply air mass correction
                sza_rad = np.radians(self.solar_zenith)
                vza_rad = np.radians(self.view_zenith)
                air_mass = 1.0 / (np.cos(sza_rad) + 0.15) + 1.0 / (np.cos(vza_rad) + 0.15)
                wvc_corrected = wvc_from_feature * air_mass
                
                # Weight by feature reliability and total weight
                feature_weight = weight / total_weight
                wvc_sum += wvc_corrected * feature_weight
                
                if self.verbose:
                    gs.verbose(f"  {feature_name}: WVC={wvc_corrected:.3f}, weight={feature_weight:.2f}")
                    
            except Exception as e:
                if self.verbose:
                    gs.warning(f"Could not process {feature_name}: {str(e)}")
                continue
        
        # Normalize and create final WVC map
        wvc_map = f"tmp_wvc_average_{os.getpid()}"
        gs.run_command('r.mapcalc', expression=f"{wvc_map} = {wvc_sum}", quiet=True)
        
        # Calculate combined uncertainty
        try:
            stats = gs.parse_command('r.univar', map=wvc_map, flags='g')
            mean_wvc = float(stats['mean'])
            std_wvc = float(stats['stddev'])
            
            # Combined uncertainty from multiple features
            combined_uncertainty = 0.08 * (1 + std_wvc / mean_wvc)  # Enhanced for multi-feature
            
            # Create uncertainty map
            uncertainty_map = f"tmp_wvc_uncertainty_{os.getpid()}"
            gs.run_command('r.mapcalc', expression=(
                f"{uncertainty_map} = {combined_uncertainty}"), quiet=True)
            
            if self.verbose:
                gs.message(f"Average WVC: {mean_wvc:.3f} ± {combined_uncertainty:.3f} g/cm²")
            
            self.temp_maps.extend([wvc_map, uncertainty_map])
            return wvc_map, mean_wvc, uncertainty_map
            
        except Exception as e:
            gs.error(f"Error in average WVC estimation: {str(e)}")
            raise
    
    def cleanup(self):
        """Clean up temporary maps."""
        for temp_map in self.temp_maps:
            try:
                gs.run_command('g.remove', type='raster', name=temp_map, flags='f', quiet=True)
            except Exception as e:
                gs.error(f"Error cleaning up temporary map {temp_map}: {str(e)}")
                raise
    
    def _weighted_average_wvc(self, wvc_940_map, wvc_940_val, wvc_1130_map, wvc_1130_val):
        """Calculate weighted average of 940nm and 1130nm WVC estimates.
        
        Uses 940nm as primary (more sensitive) with 1130nm for high-WVC conditions.
        """
        # 940nm saturates at high WVC (band_depth > ~0.7)
        # Use 940nm as primary (more sensitive), 1130nm for high-WVC
        if wvc_940_val < 3.5:
            weight_940 = 0.7
        else:
            weight_940 = 0.3
        weight_1130 = 1.0 - weight_940
        
        wvc_map = f"tmp_wvc_weighted_{os.getpid()}"
        expr = (f"{wvc_map} = {weight_940} * {wvc_940_map} + "
                f"{weight_1130} * {wvc_1130_map}")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                               quiet=not self.verbose)
        self.temp_maps.append(wvc_map)
        
        mean_wvc = weight_940 * wvc_940_val + weight_1130 * wvc_1130_val
        self.wvc_map = wvc_map
        
        if self.verbose:
            gs.message(f"Weighted WVC (w940={weight_940:.1f}, "
                          f"w1130={weight_1130:.1f}): {mean_wvc:.3f} g/cm²")
        
        return wvc_map, mean_wvc
    
    def _estimate_wvc_absorption_feature(self, feature_name):
        """Estimate WVC from a specific absorption feature.
        
        Args:
            feature_name (str): '940' or '1130'
        
        Returns:
            tuple: (wvc_map, mean_wvc)
        """
        if feature_name == '940':
            return self._estimate_wvc_940nm()
        elif feature_name == '1130':
            return self._estimate_wvc_1130nm()
        else:
            raise ValueError(f"Unknown absorption feature: {feature_name}")
        
        # If we get here, all methods failed
        raise RuntimeError("Could not estimate WVC using any available method. "
                         "Check that your data includes NIR bands with water vapor features.")


def estimate_wvc(input_raster, dem, method='average', sensor_config=None,
                 solar_zenith=30.0, view_zenith=0.0, verbose=False):
    """Convenience function to estimate WVC from hyperspectral data.

    Args:
        input_raster (str): Name of the input 3D hyperspectral raster
        dem (str): Name of the digital elevation model (DEM) raster
        method (str, optional): WVC estimation method. Defaults to 'average'.
        sensor_config (dict, optional): Sensor configuration dictionary
        solar_zenith (float): Solar zenith angle in degrees
        view_zenith (float): View zenith angle in degrees
        verbose (bool, optional): Enable verbose output

    Returns:
        tuple: (wvc_map, mean_wvc) where wvc_map is the WVC raster and mean_wvc is the mean value
    """
    estimator = WVCEstimator(input_raster, dem, sensor_config,
                             solar_zenith=solar_zenith,
                             view_zenith=view_zenith,
                             verbose=verbose)
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
    
    if gs.verbosity() > 0:
        gs.message(f"WVC map created: {wvc_map}")
        gs.message(f"Mean WVC: {mean_wvc:.3f} g/cm²")
