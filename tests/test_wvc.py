#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for WVC (Water Vapor Content) estimation algorithm.

These tests verify the core WVC math without requiring GRASS GIS.
All algorithms are reimplemented inline from lib/wvc.py.
"""

import unittest
import sys
import os
import math
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

# Constants copied from lib/wvc.py (cannot import directly due to GRASS dependency)
WATER_VAPOR_BANDS = {
    '940': {
        'center': 940,
        'window': (860, 1010),
        'absorption': (925, 955)
    },
    '1130': {
        'center': 1130,
        'window': (1060, 1200),
        'absorption': (1110, 1150)
    }
}
COEFS_940 = {'scale': 28.0, 'offset': 0.0}
COEFS_1130 = {'scale': 42.0, 'offset': 0.0}
DEFAULT_WVC_MIN = 0.1
DEFAULT_WVC_MAX = 8.0


# --- Reimplemented pure-math functions (no GRASS dependency) ---

def compute_airmass(solar_zenith, view_zenith):
    """Two-way air mass factor."""
    return (1.0 / math.cos(math.radians(solar_zenith)) +
            1.0 / math.cos(math.radians(view_zenith)))


def continuum_removal_depth(abs_val, left_val, right_val, weight):
    """Band depth via continuum removal.

    continuum = left + (right - left) * weight
    depth = 1 - abs / continuum
    """
    continuum = left_val + (right_val - left_val) * weight
    if continuum <= 0:
        return 0.0
    return 1.0 - abs_val / continuum


def interpolation_weight(wl_abs, wl_left, wl_right):
    """Linear interpolation weight for continuum."""
    return (wl_abs - wl_left) / (wl_right - wl_left)


def band_depth_to_wvc(band_depth, airmass, coefs):
    """Convert airmass-normalized band depth to WVC."""
    return (band_depth / airmass) * coefs['scale'] + coefs['offset']


def weighted_combination(wvc_940, wvc_1130, threshold=3.5):
    """Weighted combination of 940nm and 1130nm WVC estimates."""
    if wvc_940 < threshold:
        w940 = 0.7
    else:
        w940 = 0.3
    w1130 = 1.0 - w940
    return w940 * wvc_940 + w1130 * wvc_1130, w940, w1130


def select_nearest_bands(bands, target_wl, n=3):
    """Select the n bands closest to a target wavelength."""
    sorted_bands = sorted(bands, key=lambda b: abs(b['wavelength'] - target_wl))
    return sorted_bands[:min(n, len(sorted_bands))]


def clamp_wvc(wvc, wvc_min=DEFAULT_WVC_MIN, wvc_max=DEFAULT_WVC_MAX):
    """Clamp WVC to valid range."""
    return max(wvc_min, min(wvc, wvc_max))


# --- Test classes ---

class TestAirMassFactor(unittest.TestCase):
    """Air mass factor computation."""

    def test_nadir_sun_nadir_view(self):
        """SZA=0, VZA=0: airmass = 1/1 + 1/1 = 2.0."""
        am = compute_airmass(0, 0)
        self.assertAlmostEqual(am, 2.0, places=5)

    def test_sza60_nadir(self):
        """SZA=60, VZA=0: airmass = 1/cos(60) + 1 = 2 + 1 = 3.0."""
        am = compute_airmass(60, 0)
        self.assertAlmostEqual(am, 3.0, places=5)

    def test_symmetry(self):
        """Swapping SZA and VZA gives the same airmass."""
        am1 = compute_airmass(30, 15)
        am2 = compute_airmass(15, 30)
        self.assertAlmostEqual(am1, am2, places=5)

    def test_large_sza(self):
        """Large SZA (70 deg) increases airmass substantially."""
        am = compute_airmass(70, 0)
        # 1/cos(70) ~ 2.92, plus 1.0 = 3.92
        self.assertGreater(am, 3.5)
        self.assertLess(am, 4.5)

    def test_increases_with_zenith(self):
        """Airmass increases as zenith angles grow."""
        am_0 = compute_airmass(0, 0)
        am_30 = compute_airmass(30, 0)
        am_60 = compute_airmass(60, 0)
        self.assertLess(am_0, am_30)
        self.assertLess(am_30, am_60)


class TestBandDepth(unittest.TestCase):
    """Continuum-removed band depth."""

    def test_no_absorption(self):
        """No absorption (abs == continuum): depth = 0."""
        # All values equal -> continuum = abs -> depth = 0
        depth = continuum_removal_depth(1.0, 1.0, 1.0, 0.5)
        self.assertAlmostEqual(depth, 0.0, places=5)

    def test_full_absorption(self):
        """Full absorption (abs = 0): depth = 1."""
        depth = continuum_removal_depth(0.0, 1.0, 1.0, 0.5)
        self.assertAlmostEqual(depth, 1.0, places=5)

    def test_partial_absorption(self):
        """Partial absorption: 0 < depth < 1."""
        depth = continuum_removal_depth(0.7, 1.0, 1.0, 0.5)
        self.assertGreater(depth, 0)
        self.assertLess(depth, 1)
        self.assertAlmostEqual(depth, 0.3, places=5)

    def test_interpolation_weight(self):
        """Weight = (wl_abs - wl_left) / (wl_right - wl_left)."""
        w = interpolation_weight(940, 860, 1010)
        expected = (940 - 860) / (1010 - 860)
        self.assertAlmostEqual(w, expected, places=5)
        self.assertGreater(w, 0)
        self.assertLess(w, 1)

    def test_weighted_continuum(self):
        """Continuum interpolated between left and right shoulders."""
        left = 1.0
        right = 0.8
        weight = 0.5
        # continuum = 1.0 + (0.8 - 1.0) * 0.5 = 0.9
        depth = continuum_removal_depth(0.7, left, right, weight)
        continuum = left + (right - left) * weight
        self.assertAlmostEqual(continuum, 0.9, places=5)
        self.assertAlmostEqual(depth, 1.0 - 0.7 / 0.9, places=5)


class TestAirMassNormalization(unittest.TestCase):
    """Air mass normalization of band depth."""

    def test_nadir_normalization(self):
        """SZA=0, VZA=0 (airmass=2): vertical depth = depth/2."""
        airmass = compute_airmass(0, 0)
        band_depth = 0.4
        vertical_depth = band_depth / airmass
        self.assertAlmostEqual(vertical_depth, 0.2, places=5)

    def test_sza60_normalization(self):
        """SZA=60, VZA=0 (airmass=3): vertical depth = depth/3."""
        airmass = compute_airmass(60, 0)
        band_depth = 0.6
        vertical_depth = band_depth / airmass
        self.assertAlmostEqual(vertical_depth, 0.2, places=5)

    def test_higher_airmass_smaller_vertical(self):
        """Higher airmass -> dividing recovers smaller vertical depth."""
        band_depth = 0.5
        am_low = compute_airmass(0, 0)    # 2.0
        am_high = compute_airmass(60, 0)  # 3.0
        self.assertGreater(band_depth / am_low, band_depth / am_high)


class TestWVCCoefficients(unittest.TestCase):
    """WVC conversion coefficients."""

    def test_coefs_940_values(self):
        """COEFS_940 should have scale=28.0, offset=0.0."""
        self.assertEqual(COEFS_940['scale'], 28.0)
        self.assertEqual(COEFS_940['offset'], 0.0)

    def test_coefs_1130_values(self):
        """COEFS_1130 should have scale=42.0, offset=0.0."""
        self.assertEqual(COEFS_1130['scale'], 42.0)
        self.assertEqual(COEFS_1130['offset'], 0.0)

    def test_zero_band_depth(self):
        """Zero band depth -> WVC = 0."""
        wvc = band_depth_to_wvc(0.0, 2.0, COEFS_940)
        self.assertAlmostEqual(wvc, 0.0, places=5)

    def test_wvc_formula(self):
        """WVC = (band_depth / airmass) * scale."""
        band_depth = 0.2
        airmass = 2.0
        wvc = band_depth_to_wvc(band_depth, airmass, COEFS_940)
        expected = (0.2 / 2.0) * 28.0
        self.assertAlmostEqual(wvc, expected, places=5)
        self.assertAlmostEqual(wvc, 2.8, places=5)

    def test_typical_wvc_range(self):
        """Known band depth -> expected WVC in 1-5 g/cm² range."""
        # Typical observed band depth ~0.15-0.30 at airmass ~2
        for bd in [0.15, 0.20, 0.25, 0.30]:
            wvc_940 = band_depth_to_wvc(bd, 2.0, COEFS_940)
            self.assertGreater(wvc_940, 1.0)
            self.assertLess(wvc_940, 5.0)

    def test_1130_larger_scale(self):
        """1130nm has larger scale factor -> higher WVC for same depth."""
        bd = 0.2
        airmass = 2.0
        wvc_940 = band_depth_to_wvc(bd, airmass, COEFS_940)
        wvc_1130 = band_depth_to_wvc(bd, airmass, COEFS_1130)
        self.assertGreater(wvc_1130, wvc_940)


class TestWeightedCombination(unittest.TestCase):
    """Weighted combination of 940nm and 1130nm WVC estimates."""

    def test_normal_wvc_weights(self):
        """Normal WVC (<3.5): w940=0.7, w1130=0.3."""
        wvc, w940, w1130 = weighted_combination(2.0, 2.5)
        self.assertAlmostEqual(w940, 0.7, places=5)
        self.assertAlmostEqual(w1130, 0.3, places=5)

    def test_high_wvc_weights(self):
        """High WVC (>=3.5): w940=0.3, w1130=0.7."""
        wvc, w940, w1130 = weighted_combination(4.0, 3.5)
        self.assertAlmostEqual(w940, 0.3, places=5)
        self.assertAlmostEqual(w1130, 0.7, places=5)

    def test_weights_sum_to_one(self):
        """Weights always sum to 1.0."""
        for wvc_940_val in [1.0, 2.0, 3.0, 3.5, 4.0, 5.0]:
            _, w940, w1130 = weighted_combination(wvc_940_val, 2.0)
            self.assertAlmostEqual(w940 + w1130, 1.0, places=5)

    def test_result_formula(self):
        """result = w940 * wvc_940 + w1130 * wvc_1130."""
        wvc_940_val = 2.0
        wvc_1130_val = 2.5
        wvc, w940, w1130 = weighted_combination(wvc_940_val, wvc_1130_val)
        expected = w940 * wvc_940_val + w1130 * wvc_1130_val
        self.assertAlmostEqual(wvc, expected, places=5)

    def test_threshold_boundary(self):
        """At threshold=3.5 exactly, high-WVC weights apply."""
        _, w940, w1130 = weighted_combination(3.5, 3.0)
        self.assertAlmostEqual(w940, 0.3, places=5)
        self.assertAlmostEqual(w1130, 0.7, places=5)

    def test_just_below_threshold(self):
        """Just below threshold, normal weights apply."""
        _, w940, w1130 = weighted_combination(3.49, 3.0)
        self.assertAlmostEqual(w940, 0.7, places=5)
        self.assertAlmostEqual(w1130, 0.3, places=5)


class TestBandSelection(unittest.TestCase):
    """_select_nearest_bands() logic."""

    def setUp(self):
        """Create a mock band list with ~5nm spacing."""
        self.bands = [{'band': i, 'wavelength': 900 + i * 5}
                      for i in range(30)]  # 900-1045nm

    def test_selects_nearest(self):
        """Selects n closest bands to target wavelength."""
        selected = select_nearest_bands(self.bands, 940, n=3)
        self.assertEqual(len(selected), 3)
        # All selected should be near 940nm
        for b in selected:
            self.assertLess(abs(b['wavelength'] - 940), 10)

    def test_returns_min_n_available(self):
        """Returns min(n, available) bands."""
        small_list = [{'band': 1, 'wavelength': 935},
                      {'band': 2, 'wavelength': 940}]
        selected = select_nearest_bands(small_list, 940, n=5)
        self.assertEqual(len(selected), 2)

    def test_sorted_by_distance(self):
        """Result is sorted by distance to target."""
        selected = select_nearest_bands(self.bands, 940, n=5)
        distances = [abs(b['wavelength'] - 940) for b in selected]
        self.assertEqual(distances, sorted(distances))

    def test_single_band(self):
        """Selecting 1 band returns the closest."""
        selected = select_nearest_bands(self.bands, 942, n=1)
        self.assertEqual(len(selected), 1)
        # Should pick 940 (index 8: 900 + 8*5 = 940)
        self.assertEqual(selected[0]['wavelength'], 940)

    def test_empty_input(self):
        """Empty band list returns empty result."""
        selected = select_nearest_bands([], 940, n=3)
        self.assertEqual(len(selected), 0)


class TestAbsorptionWindows(unittest.TestCase):
    """Absorption window definitions."""

    def test_940nm_window(self):
        """940nm window: (860, 1010), absorption: (925, 955)."""
        feature = WATER_VAPOR_BANDS['940']
        self.assertEqual(feature['window'], (860, 1010))
        self.assertEqual(feature['absorption'], (925, 955))
        self.assertEqual(feature['center'], 940)

    def test_1130nm_window(self):
        """1130nm window: (1060, 1200), absorption: (1110, 1150)."""
        feature = WATER_VAPOR_BANDS['1130']
        self.assertEqual(feature['window'], (1060, 1200))
        self.assertEqual(feature['absorption'], (1110, 1150))
        self.assertEqual(feature['center'], 1130)

    def test_left_shoulder_before_absorption(self):
        """Left shoulder bands have wavelength < absorption min."""
        for name, feature in WATER_VAPOR_BANDS.items():
            window_min = feature['window'][0]
            abs_min = feature['absorption'][0]
            self.assertLess(window_min, abs_min,
                            msg=f"{name}nm: window_min >= abs_min")

    def test_right_shoulder_after_absorption(self):
        """Right shoulder bands have wavelength > absorption max."""
        for name, feature in WATER_VAPOR_BANDS.items():
            window_max = feature['window'][1]
            abs_max = feature['absorption'][1]
            self.assertGreater(window_max, abs_max,
                               msg=f"{name}nm: window_max <= abs_max")

    def test_tanager_band_coverage(self):
        """With ~5nm spacing, many bands should fall in each region."""
        for name, feature in WATER_VAPOR_BANDS.items():
            window_min, window_max = feature['window']
            abs_min, abs_max = feature['absorption']

            # Simulate Tanager bands (~5nm spacing)
            tanager_bands = [{'wavelength': wl, 'band': i}
                             for i, wl in enumerate(
                                 np.arange(376, 2500, 5.5))]

            left = [b for b in tanager_bands
                    if window_min <= b['wavelength'] < abs_min]
            absorb = [b for b in tanager_bands
                      if abs_min <= b['wavelength'] <= abs_max]
            right = [b for b in tanager_bands
                     if abs_max < b['wavelength'] <= window_max]

            self.assertGreater(len(left), 3,
                               msg=f"{name}nm left shoulder: {len(left)} bands")
            self.assertGreater(len(absorb), 2,
                               msg=f"{name}nm absorption: {len(absorb)} bands")
            self.assertGreater(len(right), 3,
                               msg=f"{name}nm right shoulder: {len(right)} bands")


class TestWVCBounds(unittest.TestCase):
    """WVC clamping to valid range [0.1, 8.0]."""

    def test_normal_value(self):
        """Normal WVC passes through."""
        self.assertEqual(clamp_wvc(2.5), 2.5)
        self.assertEqual(clamp_wvc(5.0), 5.0)

    def test_below_min(self):
        """WVC below minimum is clamped to 0.1."""
        self.assertEqual(clamp_wvc(0.0), DEFAULT_WVC_MIN)
        self.assertEqual(clamp_wvc(-1.0), DEFAULT_WVC_MIN)

    def test_above_max(self):
        """WVC above maximum is clamped to 8.0."""
        self.assertEqual(clamp_wvc(10.0), DEFAULT_WVC_MAX)
        self.assertEqual(clamp_wvc(100.0), DEFAULT_WVC_MAX)

    def test_boundary_values(self):
        """Boundary values are preserved."""
        self.assertEqual(clamp_wvc(DEFAULT_WVC_MIN), DEFAULT_WVC_MIN)
        self.assertEqual(clamp_wvc(DEFAULT_WVC_MAX), DEFAULT_WVC_MAX)


if __name__ == '__main__':
    unittest.main()
