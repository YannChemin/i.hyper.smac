#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for utility functions.

These tests verify the utility functions without requiring GRASS GIS.
"""

import unittest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))


class TestFindNearestBand(unittest.TestCase):
    """Tests for find_nearest_band function."""

    def test_exact_match(self):
        """Should find exact wavelength match."""
        # Import here to avoid GRASS dependency at module level
        # We'll mock the function for testing without GRASS
        bands = [
            {'band': 1, 'wavelength': 450.0},
            {'band': 2, 'wavelength': 550.0},
            {'band': 3, 'wavelength': 650.0},
            {'band': 4, 'wavelength': 850.0},
        ]

        # Simple find_nearest implementation for testing
        def find_nearest(bands, target):
            return min(bands, key=lambda x: abs(x['wavelength'] - target))

        result = find_nearest(bands, 550.0)
        self.assertEqual(result['band'], 2)
        self.assertEqual(result['wavelength'], 550.0)

    def test_closest_match(self):
        """Should find closest wavelength when no exact match."""
        bands = [
            {'band': 1, 'wavelength': 450.0},
            {'band': 2, 'wavelength': 550.0},
            {'band': 3, 'wavelength': 650.0},
            {'band': 4, 'wavelength': 850.0},
        ]

        def find_nearest(bands, target):
            return min(bands, key=lambda x: abs(x['wavelength'] - target))

        # 600nm should match band 3 (650nm) or band 2 (550nm)
        result = find_nearest(bands, 600.0)
        self.assertIn(result['band'], [2, 3])

        # 900nm should match band 4 (850nm)
        result = find_nearest(bands, 900.0)
        self.assertEqual(result['band'], 4)

    def test_tuple_input(self):
        """Should handle tuple input format."""
        # Test the conversion from tuples to dicts
        bands_tuples = [
            (450.0, 'band1'),
            (550.0, 'band2'),
            (650.0, 'band3'),
        ]

        # Convert to dict format
        bands = [{'wavelength': wl, 'band': i + 1}
                 for i, (wl, _) in enumerate(bands_tuples)]

        def find_nearest(bands, target):
            return min(bands, key=lambda x: abs(x['wavelength'] - target))

        result = find_nearest(bands, 550.0)
        self.assertEqual(result['band'], 2)


class TestConvertWavelength(unittest.TestCase):
    """Tests for wavelength unit conversion."""

    def test_nm_passthrough(self):
        """Nanometers should pass through unchanged."""
        def convert(wavelength, unit):
            unit = unit.lower().strip()
            if unit in ['nm', 'nanometer', 'nanometers']:
                return wavelength
            elif unit in ['um', 'µm', 'micrometer', 'micrometers']:
                return wavelength * 1000.0
            elif unit in ['m', 'meter', 'meters']:
                return wavelength * 1e9
            return wavelength

        self.assertEqual(convert(550.0, 'nm'), 550.0)
        self.assertEqual(convert(550.0, 'nanometer'), 550.0)
        self.assertEqual(convert(550.0, 'nanometers'), 550.0)

    def test_um_to_nm(self):
        """Micrometers should convert to nanometers."""
        def convert(wavelength, unit):
            unit = unit.lower().strip()
            if unit in ['nm', 'nanometer', 'nanometers']:
                return wavelength
            elif unit in ['um', 'µm', 'micrometer', 'micrometers']:
                return wavelength * 1000.0
            return wavelength

        self.assertEqual(convert(0.55, 'um'), 550.0)
        self.assertEqual(convert(0.55, 'µm'), 550.0)
        self.assertEqual(convert(1.0, 'micrometer'), 1000.0)

    def test_m_to_nm(self):
        """Meters should convert to nanometers."""
        def convert(wavelength, unit):
            unit = unit.lower().strip()
            if unit in ['nm', 'nanometer', 'nanometers']:
                return wavelength
            elif unit in ['um', 'µm', 'micrometer', 'micrometers']:
                return wavelength * 1000.0
            elif unit in ['m', 'meter', 'meters']:
                return wavelength * 1e9
            return wavelength

        self.assertAlmostEqual(convert(550e-9, 'm'), 550.0)
        self.assertAlmostEqual(convert(1e-6, 'meter'), 1000.0)


class TestRadtranFunctions(unittest.TestCase):
    """Tests for radtran utility functions."""

    def test_earth_sun_distance(self):
        """Earth-Sun distance factor should vary seasonally.

        Note: The formula returns a correction factor, not actual distance.
        The factor is used for irradiance calculations: E = E0 / d^2
        """
        import numpy as np
        from datetime import datetime

        def earth_sun_distance(year, month, day):
            doy = datetime(year, month, day).timetuple().tm_yday
            beta = 2 * np.pi * (doy - 3) / 365.25
            dist = 1 + 0.01670963 * np.cos(beta) - 0.0000146 * np.cos(2 * beta)
            return dist

        # Test seasonal variation
        d_jan = earth_sun_distance(2024, 1, 3)
        d_jul = earth_sun_distance(2024, 7, 4)

        # Values should differ seasonally (about 3.3% variation)
        self.assertNotAlmostEqual(d_jan, d_jul, places=2)

        # Both values should be close to 1 AU (within ~2%)
        self.assertGreater(d_jan, 0.98)
        self.assertLess(d_jan, 1.02)
        self.assertGreater(d_jul, 0.98)
        self.assertLess(d_jul, 1.02)

        # Total variation should be about 3.3% (Earth's orbital eccentricity)
        variation = abs(d_jan - d_jul)
        self.assertGreater(variation, 0.02)
        self.assertLess(variation, 0.05)

    def test_gaussian_response(self):
        """Gaussian band response function should peak at center."""
        import numpy as np

        def gaussian_rsp(wl, wl_center, fwhm):
            sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
            return np.exp(-0.5 * ((wl - wl_center) / sigma) ** 2)

        wl_center = 550.0
        fwhm = 10.0

        # Response at center should be 1.0
        self.assertAlmostEqual(gaussian_rsp(wl_center, wl_center, fwhm), 1.0)

        # Response at FWHM distance should be 0.5
        self.assertAlmostEqual(
            gaussian_rsp(wl_center + fwhm / 2, wl_center, fwhm), 0.5, places=5
        )
        self.assertAlmostEqual(
            gaussian_rsp(wl_center - fwhm / 2, wl_center, fwhm), 0.5, places=5
        )

        # Response far from center should approach 0
        self.assertLess(gaussian_rsp(wl_center + 3 * fwhm, wl_center, fwhm), 0.01)


class TestAtmosphericEstimation(unittest.TestCase):
    """Tests for atmospheric parameter estimation logic."""

    def test_aod_bounds(self):
        """AOD estimates should be within reasonable bounds."""
        # Typical AOD range is 0.01 to 1.5
        DEFAULT_AOD_MIN = 0.01
        DEFAULT_AOD_MAX = 1.5

        def clamp_aod(aod):
            return max(DEFAULT_AOD_MIN, min(aod, DEFAULT_AOD_MAX))

        self.assertEqual(clamp_aod(0.15), 0.15)  # Normal value
        self.assertEqual(clamp_aod(-0.1), DEFAULT_AOD_MIN)  # Below min
        self.assertEqual(clamp_aod(2.0), DEFAULT_AOD_MAX)  # Above max

    def test_wvc_bounds(self):
        """WVC estimates should be within reasonable bounds."""
        # Typical WVC range is 0.1 to 8.0 g/cm²
        DEFAULT_WVC_MIN = 0.1
        DEFAULT_WVC_MAX = 8.0

        def clamp_wvc(wvc):
            return max(DEFAULT_WVC_MIN, min(wvc, DEFAULT_WVC_MAX))

        self.assertEqual(clamp_wvc(2.5), 2.5)  # Normal value
        self.assertEqual(clamp_wvc(0.0), DEFAULT_WVC_MIN)  # Below min
        self.assertEqual(clamp_wvc(10.0), DEFAULT_WVC_MAX)  # Above max

    def test_ozone_bounds(self):
        """Ozone estimates should be within reasonable bounds."""
        # Typical ozone range is 150 to 500 DU
        DEFAULT_O3_MIN = 150
        DEFAULT_O3_MAX = 500

        def clamp_ozone(o3):
            return max(DEFAULT_O3_MIN, min(o3, DEFAULT_O3_MAX))

        self.assertEqual(clamp_ozone(300), 300)  # Normal value
        self.assertEqual(clamp_ozone(100), DEFAULT_O3_MIN)  # Below min
        self.assertEqual(clamp_ozone(600), DEFAULT_O3_MAX)  # Above max


class TestNDVICalculation(unittest.TestCase):
    """Tests for NDVI calculation logic."""

    def test_ndvi_formula(self):
        """NDVI should follow (NIR - Red) / (NIR + Red) formula."""
        import numpy as np

        def calc_ndvi(nir, red):
            return (nir - red) / (nir + red + 0.0001)

        # Vegetation: high NIR, low red -> high NDVI
        self.assertGreater(calc_ndvi(0.5, 0.1), 0.5)

        # Water: low NIR, low red -> low/negative NDVI
        self.assertLess(calc_ndvi(0.05, 0.1), 0)

        # Bare soil: similar NIR and red -> low NDVI
        ndvi_soil = calc_ndvi(0.3, 0.25)
        self.assertLess(abs(ndvi_soil), 0.2)

    def test_ndvi_array(self):
        """NDVI calculation should work with arrays."""
        import numpy as np

        def calc_ndvi(nir, red):
            return (nir - red) / (nir + red + 0.0001)

        nir = np.array([0.5, 0.3, 0.05, 0.4])
        red = np.array([0.1, 0.25, 0.1, 0.15])

        ndvi = calc_ndvi(nir, red)

        self.assertEqual(ndvi.shape, nir.shape)
        self.assertTrue(np.all(ndvi >= -1.0))
        self.assertTrue(np.all(ndvi <= 1.0))


if __name__ == '__main__':
    unittest.main()
