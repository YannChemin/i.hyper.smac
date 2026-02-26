#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AOD (Aerosol Optical Depth) estimation module for hyperspectral imagery.

Implements the MODIS Dark Target algorithm (Kaufman 1997, Levy 2007) adapted
for hyperspectral data, using VIS/SWIR surface reflectance relationships
with single-scattering AOD inversion.

Fallback chain: Dark Target -> Dark Pixel -> Constant(0.15)
"""

import os
import math
import grass.script as gs

try:
    from .utils import get_band_info, find_nearest_band
except ImportError:
    from utils import get_band_info, find_nearest_band

try:
    from . import radtran
except ImportError:
    import radtran

# Target wavelengths for Dark Target algorithm (nm)
BLUE_WL = 470.0
RED_WL = 660.0
NIR_WL = 860.0
SWIR_WL = 2130.0

# Continental aerosol parameters
OMEGA_0 = 0.89        # Single scattering albedo
ASYM_G = 0.65         # Asymmetry parameter (Henyey-Greenstein)
ANGSTROM_DEFAULT = 1.3 # Default Angstrom exponent


class AODEstimator:
    """Estimate Aerosol Optical Depth using the Dark Target algorithm."""

    def __init__(self, input_raster, dem, solar_zenith, view_zenith,
                 solar_azimuth, view_azimuth, pressure=1013.25,
                 verbose=False):
        self.input_raster = input_raster
        self.dem = dem
        self.solar_zenith = solar_zenith
        self.view_zenith = view_zenith
        self.solar_azimuth = solar_azimuth
        self.view_azimuth = view_azimuth
        self.pressure = pressure
        self.verbose = verbose
        self.band_info = []
        self.temp_maps = []
        self.retrieved_alpha = None  # Angstrom exponent from Dark Target

        # Precompute cosines (angles in degrees -> radians)
        self.mu_s = math.cos(math.radians(solar_zenith))
        self.mu_v = math.cos(math.radians(view_zenith))

        # Scattering angle
        theta_s_r = math.radians(solar_zenith)
        theta_v_r = math.radians(view_zenith)
        dphi_r = math.radians(solar_azimuth - view_azimuth)
        cos_scatter = (-math.cos(theta_s_r) * math.cos(theta_v_r)
                       + math.sin(theta_s_r) * math.sin(theta_v_r)
                       * math.cos(dphi_r))
        cos_scatter = max(-1.0, min(1.0, cos_scatter))
        self.scatter_angle = math.acos(cos_scatter)

        # Henyey-Greenstein phase function
        g = ASYM_G
        cos_scat = math.cos(self.scatter_angle)
        self.phase_hg = ((1 - g**2)
                         / (1 + g**2 - 2 * g * cos_scat)**1.5)

        self._get_band_info()

    def _get_band_info(self):
        self.band_info = get_band_info(self.input_raster, verbose=self.verbose)

    def _find_nearest_band(self, target_wavelength):
        return find_nearest_band(self.band_info, target_wavelength)

    def extract_band_to_2d(self, band_num, output_map=None):
        """Extract a single band from 3D raster to 2D raster."""
        import time
        import random

        if not output_map:
            output_map = (f"tmp_aod_band_{band_num}_{int(time.time())}"
                          f"_{random.randint(1000, 9999)}")

        try:
            gs.run_command('g.remove', flags='f', type='raster',
                           pattern=f"{output_map}*", quiet=True)
            gs.run_command('g.region', t=band_num + 0.1, b=band_num,
                           quiet=True)
            gs.run_command('r3.to.rast', input=self.input_raster,
                           output=output_map, overwrite=True,
                           quiet=not self.verbose)

            output_file = f"{output_map}_00001"
            gs.run_command('g.rename', raster=f"{output_file},{output_map}",
                           overwrite=True, quiet=True)
            gs.run_command('g.region', raster_3d=self.input_raster,
                           quiet=True)

            if self.verbose:
                gs.message(f"Extracted band {band_num} to {output_map}")

            self.temp_maps.append(output_map)
            self.temp_maps.append(output_file)
            return output_map

        except Exception as e:
            gs.warning(f"Error extracting band {band_num}: {e}")
            raise
        finally:
            gs.run_command('g.remove', flags='f', type='raster',
                           pattern=f"{output_map}_00001", quiet=True)
            gs.run_command('g.remove', flags='f', type='raster_3d',
                           name='RASTER3D_MASK', quiet=True)

    def _compute_E0(self, wavelength, fwhm):
        """Compute exo-atmospheric irradiance with analytical fallback.

        Returns E0 in mW/(m^2 nm).
        """
        try:
            e0 = radtran.E0(wavelength, fwhm, verbose=False)
            if e0 is not None and e0 > 0:
                return e0
        except Exception:
            pass

        # Analytical blackbody fallback (Planck function at 5778K)
        wl_m = wavelength * 1e-9
        h = 6.62607015e-34
        c = 2.998e8
        k = 1.381e-23
        T = 5778.0
        hc_kT = h * c / (wl_m * k * T)
        B = 2 * h * c**2 / (wl_m**5 * (math.exp(hc_kT) - 1.0))
        # Solid angle of sun as seen from 1 AU: Omega_sun = pi * (R_sun/d)^2
        # R_sun/d_AU ~ 0.00465 rad => Omega = pi * 0.00465^2 = 6.794e-5 sr
        # E0 = Omega * B, convert W/(m^2 m) to mW/(m^2 nm): * 1e-6
        e0 = 6.794e-5 * B * 1e-6
        if self.verbose:
            gs.message(f"  E0({wavelength:.0f}nm) = {e0:.4f} mW/(m^2 nm) "
                       f"[analytical fallback]")
        return e0

    def _radiance_to_reflectance(self, radiance_map, wavelength, fwhm):
        """Convert radiance raster to TOA reflectance.

        Input radiance: W/(m^2 sr um) = mW/(m^2 sr nm)
        E0: mW/(m^2 nm)
        Output: unitless TOA reflectance
        """
        e0 = self._compute_E0(wavelength, fwhm)

        # Earth-Sun distance (assume mid-year if no timestamp)
        try:
            info = gs.raster3d_info(self.input_raster)
            ts = info.get('timestamp', '')
            if ts and ts.strip():
                from datetime import datetime
                dt = datetime.strptime(ts.strip()[:10], '%Y-%m-%d')
                d = radtran.earth_sun_distance(dt.year, dt.month, dt.day)
            else:
                d = 1.0
        except Exception:
            d = 1.0

        # rho = pi * L * d^2 / (E0 * cos(theta_s))
        scale = math.pi * d**2 / (e0 * self.mu_s)

        refl_map = f"tmp_aod_refl_{int(wavelength)}_{os.getpid()}"
        expr = f"{refl_map} = {radiance_map} * {scale}"
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(refl_map)

        if self.verbose:
            gs.message(f"  Converted {radiance_map} to reflectance "
                       f"(E0={e0:.2f}, d={d:.4f}, scale={scale:.6f})")

        return refl_map

    def _estimate_aod_dark_target(self):
        """Estimate AOD using the MODIS Dark Target algorithm.

        Steps:
        1. Extract and convert 4 key bands to TOA reflectance
        2. Select Dark Target pixels using SWIR and NDVI
        3. Predict surface VIS reflectance from SWIR (Kaufman 1997)
        4. Compute atmospheric path radiance
        5. Invert to AOD using single-scattering approximation
        6. Derive Angstrom exponent from blue/red ratio
        """
        if self.verbose:
            gs.message("Dark Target AOD estimation")
            gs.message(f"  Geometry: SZA={self.solar_zenith:.1f}, "
                       f"VZA={self.view_zenith:.1f}, "
                       f"scatter={math.degrees(self.scatter_angle):.1f} deg")
            gs.message(f"  Phase function P(Theta)={self.phase_hg:.3f}")

        # Step 1: Extract and convert key bands to TOA reflectance
        bands = {}
        for name, target_wl in [('blue', BLUE_WL), ('red', RED_WL),
                                 ('nir', NIR_WL), ('swir', SWIR_WL)]:
            band = self._find_nearest_band(target_wl)
            if self.verbose:
                gs.message(f"  {name}: band {band['band']} "
                           f"({band['wavelength']:.1f} nm)")
            radiance = self.extract_band_to_2d(
                band['band'], output_map=f"tmp_aod_{name}_rad")
            refl = self._radiance_to_reflectance(
                radiance, band['wavelength'], band.get('fwhm', 10.0))
            bands[name] = {'map': refl, 'band': band}

        # Step 2: Dark Target pixel selection (Levy 2007)
        # SWIR range: 0.01 < rho_SWIR < 0.25
        swir_map = bands['swir']['map']
        red_map = bands['red']['map']
        nir_map = bands['nir']['map']
        blue_map = bands['blue']['map']

        # Compute NDVI
        ndvi = f"tmp_aod_ndvi_{os.getpid()}"
        expr = (f"{ndvi} = float({nir_map} - {red_map}) "
                f"/ float({nir_map} + {red_map} + 0.0001)")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(ndvi)

        # Initial dark target mask: SWIR range + NDVI > 0.1
        dt_mask_initial = f"tmp_aod_dt_mask_init_{os.getpid()}"
        expr = (f"{dt_mask_initial} = if({swir_map} > 0.01 && "
                f"{swir_map} < 0.25 && {ndvi} > 0.1, 1, null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(dt_mask_initial)

        # Count initial DT pixels
        try:
            stats = gs.parse_command('r.univar', map=dt_mask_initial,
                                     flags='g')
            dt_count = int(float(stats['n']))
        except Exception:
            dt_count = 0

        if self.verbose:
            gs.message(f"  Initial Dark Target pixels: {dt_count}")

        if dt_count < 100:
            if self.verbose:
                gs.message("  Insufficient DT pixels, falling back to "
                           "dark pixel method")
            return self._estimate_aod_dark_pixel(bands)

        # Filter to 20th-50th percentile of SWIR reflectance within DT mask
        swir_masked = f"tmp_aod_swir_masked_{os.getpid()}"
        gs.run_command('r.mapcalc',
                       expression=(f"{swir_masked} = if({dt_mask_initial}, "
                                   f"{swir_map}, null())"),
                       overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(swir_masked)

        try:
            quantiles = gs.read_command('r.quantile', input=swir_masked,
                                        percentiles='20,50', quiet=True)
            lines = quantiles.strip().split('\n')
            p20 = float(lines[0].split(':')[2])
            p50 = float(lines[1].split(':')[2])
        except Exception:
            # Fallback: use the full range
            p20 = 0.01
            p50 = 0.25

        dt_mask = f"tmp_aod_dt_mask_{os.getpid()}"
        expr = (f"{dt_mask} = if({dt_mask_initial} && "
                f"{swir_map} >= {p20} && {swir_map} <= {p50}, "
                f"1, null())")
        gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                       quiet=not self.verbose)
        self.temp_maps.append(dt_mask)

        try:
            stats = gs.parse_command('r.univar', map=dt_mask, flags='g')
            dt_count = int(float(stats['n']))
        except Exception:
            dt_count = 0

        if self.verbose:
            gs.message(f"  Filtered DT pixels (p20-p50): {dt_count}")

        if dt_count < 50:
            if self.verbose:
                gs.message("  Insufficient filtered DT pixels, "
                           "using initial mask")
            dt_mask = dt_mask_initial

        # Step 3: Predict surface VIS reflectance from SWIR (Kaufman 1997)
        # rho_surf(0.47um) = 0.25 * rho_SWIR(2.1um)
        # rho_surf(0.66um) = 0.50 * rho_SWIR(2.1um)
        surf_blue = f"tmp_aod_surf_blue_{os.getpid()}"
        surf_red = f"tmp_aod_surf_red_{os.getpid()}"
        gs.run_command('r.mapcalc',
                       expression=f"{surf_blue} = {swir_map} * 0.25",
                       overwrite=True, quiet=not self.verbose)
        gs.run_command('r.mapcalc',
                       expression=f"{surf_red} = {swir_map} * 0.50",
                       overwrite=True, quiet=not self.verbose)
        self.temp_maps.extend([surf_blue, surf_red])

        # Step 4: Compute atmospheric path radiance
        # rho_path = rho_TOA - rho_surf (first-order, ignoring transmission)
        path_blue = f"tmp_aod_path_blue_{os.getpid()}"
        path_red = f"tmp_aod_path_red_{os.getpid()}"
        gs.run_command(
            'r.mapcalc',
            expression=(f"{path_blue} = if({dt_mask}, "
                        f"max({blue_map} - {surf_blue}, 0.0), null())"),
            overwrite=True, quiet=not self.verbose)
        gs.run_command(
            'r.mapcalc',
            expression=(f"{path_red} = if({dt_mask}, "
                        f"max({red_map} - {surf_red}, 0.0), null())"),
            overwrite=True, quiet=not self.verbose)
        self.temp_maps.extend([path_blue, path_red])

        # Step 5: Single-scattering AOD inversion
        # tau = rho_path * 4 * mu_s * mu_v / (omega_0 * P(Theta))
        inversion_factor = 4.0 * self.mu_s * self.mu_v / (OMEGA_0 * self.phase_hg)

        tau_blue_map = f"tmp_aod_tau_blue_{os.getpid()}"
        tau_red_map = f"tmp_aod_tau_red_{os.getpid()}"
        gs.run_command(
            'r.mapcalc',
            expression=(f"{tau_blue_map} = if({dt_mask}, "
                        f"{path_blue} * {inversion_factor}, null())"),
            overwrite=True, quiet=not self.verbose)
        gs.run_command(
            'r.mapcalc',
            expression=(f"{tau_red_map} = if({dt_mask}, "
                        f"{path_red} * {inversion_factor}, null())"),
            overwrite=True, quiet=not self.verbose)
        self.temp_maps.extend([tau_blue_map, tau_red_map])

        # Step 6: Multi-wavelength Angstrom exponent
        # alpha = -ln(tau_blue/tau_red) / ln(470/660)
        try:
            stats_blue = gs.parse_command('r.univar', map=tau_blue_map,
                                          flags='g')
            stats_red = gs.parse_command('r.univar', map=tau_red_map,
                                         flags='g')
            mean_tau_blue = float(stats_blue['mean'])
            mean_tau_red = float(stats_red['mean'])

            if mean_tau_blue > 0 and mean_tau_red > 0:
                alpha = (-math.log(mean_tau_blue / mean_tau_red)
                         / math.log(BLUE_WL / RED_WL))
                if self.verbose:
                    gs.message(f"  Angstrom exponent: {alpha:.2f} "
                               f"(tau_blue={mean_tau_blue:.4f}, "
                               f"tau_red={mean_tau_red:.4f})")
                if not (0.5 <= alpha <= 2.5):
                    if self.verbose:
                        gs.message(f"  Angstrom exponent {alpha:.2f} outside "
                                   f"[0.5, 2.5], using default {ANGSTROM_DEFAULT}")
                    alpha = ANGSTROM_DEFAULT
            else:
                alpha = ANGSTROM_DEFAULT
        except Exception:
            alpha = ANGSTROM_DEFAULT

        # Store retrieved Angstrom exponent
        self.retrieved_alpha = alpha

        # Scale blue-band AOD to 550nm
        # tau_550 = tau_blue * (550/470)^(-alpha)
        scale_to_550 = (550.0 / BLUE_WL) ** (-alpha)

        aod_map = f"tmp_aod_550_{os.getpid()}"
        gs.run_command(
            'r.mapcalc',
            expression=(f"{aod_map} = if({dt_mask}, "
                        f"min(max({tau_blue_map} * {scale_to_550}, 0.01), 2.0), "
                        f"null())"),
            overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(aod_map)

        # Statistics
        try:
            stats = gs.parse_command('r.univar', map=aod_map, flags='g')
            mean_aod = float(stats['mean'])
            mean_aod = max(0.01, min(2.0, mean_aod))
            if self.verbose:
                gs.message(f"  Dark Target AOD @ 550nm: {mean_aod:.3f} "
                           f"(n={stats['n']}, "
                           f"stddev={float(stats.get('stddev', 0)):.3f}, "
                           f"alpha={alpha:.2f})")
            return aod_map, mean_aod
        except Exception as e:
            gs.warning(f"Error computing DT AOD statistics: {e}")
            return self._estimate_aod_dark_pixel(bands)

    def _estimate_aod_dark_pixel(self, bands=None):
        """Fallback: estimate AOD from dark blue pixels.

        Assumes dark blue pixels (rho_blue < 0.05) have surface rho ~ 0.01.
        Uses the same single-scattering inversion.
        """
        if self.verbose:
            gs.message("Dark Pixel fallback AOD estimation")

        try:
            if bands and 'blue' in bands:
                blue_map = bands['blue']['map']
            else:
                band = self._find_nearest_band(BLUE_WL)
                radiance = self.extract_band_to_2d(
                    band['band'], output_map='tmp_aod_dp_blue_rad')
                blue_map = self._radiance_to_reflectance(
                    radiance, band['wavelength'], band.get('fwhm', 10.0))

            # Mask: dark blue pixels with rho < 0.05
            dp_mask = f"tmp_aod_dp_mask_{os.getpid()}"
            expr = (f"{dp_mask} = if({blue_map} > 0.001 && "
                    f"{blue_map} < 0.05, 1, null())")
            gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                           quiet=not self.verbose)
            self.temp_maps.append(dp_mask)

            try:
                stats = gs.parse_command('r.univar', map=dp_mask, flags='g')
                dp_count = int(float(stats['n']))
            except Exception:
                dp_count = 0

            if self.verbose:
                gs.message(f"  Dark pixels found: {dp_count}")

            if dp_count < 50:
                if self.verbose:
                    gs.message("  Insufficient dark pixels, "
                               "using constant AOD=0.15")
                return self._create_constant_aod(0.15), 0.15

            # Path radiance: rho_path = rho_TOA_blue - 0.01
            path_blue = f"tmp_aod_dp_path_{os.getpid()}"
            expr = (f"{path_blue} = if({dp_mask}, "
                    f"max({blue_map} - 0.01, 0.0), null())")
            gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                           quiet=not self.verbose)
            self.temp_maps.append(path_blue)

            # Single-scattering inversion
            inv_factor = (4.0 * self.mu_s * self.mu_v
                          / (OMEGA_0 * self.phase_hg))
            scale_to_550 = (550.0 / BLUE_WL) ** (-ANGSTROM_DEFAULT)

            aod_map = f"tmp_aod_dp_550_{os.getpid()}"
            expr = (f"{aod_map} = if({dp_mask}, "
                    f"min(max({path_blue} * {inv_factor} * {scale_to_550}, "
                    f"0.01), 2.0), null())")
            gs.run_command('r.mapcalc', expression=expr, overwrite=True,
                           quiet=not self.verbose)
            self.temp_maps.append(aod_map)

            stats = gs.parse_command('r.univar', map=aod_map, flags='g')
            mean_aod = max(0.01, min(2.0, float(stats['mean'])))
            if self.verbose:
                gs.message(f"  Dark Pixel AOD @ 550nm: {mean_aod:.3f}")
            return aod_map, mean_aod

        except Exception as e:
            gs.warning(f"Dark pixel method failed: {e}")
            return self._create_constant_aod(0.15), 0.15

    def _create_constant_aod(self, value):
        """Create a constant AOD map."""
        aod_map = f"tmp_aod_const_{os.getpid()}"
        gs.run_command('r.mapcalc',
                       expression=f"{aod_map} = {value}",
                       overwrite=True, quiet=not self.verbose)
        self.temp_maps.append(aod_map)
        return aod_map

    def estimate_aod(self, method='auto'):
        """Estimate AOD using the specified method.

        Fallback chain: Dark Target -> Dark Pixel -> Constant(0.15)

        Args:
            method: 'auto', 'dark_target', 'dark_pixel', or 'constant'

        Returns:
            tuple: (aod_map_name, mean_aod_value)
        """
        if method == 'auto' or method == 'dark_target':
            try:
                result = self._estimate_aod_dark_target()
                if result[1] > 0:
                    return result
            except Exception as e:
                if self.verbose:
                    gs.message(f"Dark Target failed: {e}")

            if method == 'dark_target':
                gs.warning("Dark Target failed, falling back")

            try:
                return self._estimate_aod_dark_pixel()
            except Exception as e:
                if self.verbose:
                    gs.message(f"Dark Pixel failed: {e}")
                return self._create_constant_aod(0.15), 0.15

        elif method == 'dark_pixel':
            return self._estimate_aod_dark_pixel()

        elif method == 'constant':
            return self._create_constant_aod(0.15), 0.15

        else:
            raise ValueError(f"Unknown AOD method: {method}")

    def cleanup(self):
        """Clean up temporary maps."""
        for map_name in self.temp_maps:
            if gs.find_file(map_name, element='cell')['file']:
                gs.run_command('g.remove', flags='f', type='raster',
                               name=map_name, quiet=True)


def estimate_aod(input_raster, dem, method='auto',
                 solar_zenith=30.0, view_zenith=0.0,
                 solar_azimuth=0.0, view_azimuth=0.0,
                 pressure=1013.25, verbose=False):
    """Estimate AOD from hyperspectral data.

    Args:
        input_raster: Name of the input 3D hyperspectral raster
        dem: Name of the DEM raster
        method: AOD estimation method ('auto', 'dark_target', 'dark_pixel',
                'constant')
        solar_zenith: Solar zenith angle in degrees
        view_zenith: View zenith angle in degrees
        solar_azimuth: Solar azimuth angle in degrees
        view_azimuth: View azimuth angle in degrees
        pressure: Atmospheric pressure in hPa
        verbose: Enable verbose output

    Returns:
        tuple: (aod_map, mean_aod, retrieved_alpha)
            retrieved_alpha is the Angstrom exponent from Dark Target retrieval,
            or None if not available (fallback methods).
    """
    estimator = AODEstimator(
        input_raster, dem,
        solar_zenith=solar_zenith, view_zenith=view_zenith,
        solar_azimuth=solar_azimuth, view_azimuth=view_azimuth,
        pressure=pressure, verbose=verbose
    )
    try:
        aod_map, mean_aod = estimator.estimate_aod(method)

        permanent_aod = f"{input_raster}_aod_550"
        gs.run_command('g.copy', raster=f"{aod_map},{permanent_aod}",
                       overwrite=True, quiet=not verbose)

        return permanent_aod, mean_aod, estimator.retrieved_alpha
    finally:
        estimator.cleanup()
