#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for the atmospheric LUT module (lib/lut.py).

These tests verify the three-albedo extraction formulas, LUT interpolation,
and save/load round-trips without requiring libRadtran to be installed.
"""

import unittest
import os
import sys
import tempfile
import numpy as np

# Add parent lib directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from lut import extract_three_albedo, AtmosphericLUT


class TestExtractThreeAlbedo(unittest.TestCase):
    """Tests for the three-albedo extraction formulas."""

    def test_known_analytic_inputs(self):
        """Verify extraction with known R_atm, T_scat, s values.

        Forward model:  r(rho) = R_atm + T_scat * rho / (1 - s * rho)

        Given known (R_atm, T_scat, s), compute r at rho=0, 0.1, 0.5
        then verify extraction recovers the originals.
        """
        R_atm_true = 0.05
        T_scat_true = 0.80
        s_true = 0.10

        rho1, rho2 = 0.1, 0.5

        r0 = R_atm_true
        r1 = R_atm_true + T_scat_true * rho1 / (1 - s_true * rho1)
        r2 = R_atm_true + T_scat_true * rho2 / (1 - s_true * rho2)

        R_atm, T_scat, s_alb = extract_three_albedo(r0, r1, r2, rho1, rho2)

        self.assertAlmostEqual(R_atm, R_atm_true, places=10)
        self.assertAlmostEqual(T_scat, T_scat_true, places=8)
        self.assertAlmostEqual(s_alb, s_true, places=8)

    def test_zero_spherical_albedo(self):
        """When s=0, r(rho) = R_atm + T_scat * rho (linear)."""
        R_atm_true = 0.03
        T_scat_true = 0.90
        s_true = 0.0

        rho1, rho2 = 0.1, 0.5

        r0 = R_atm_true
        r1 = R_atm_true + T_scat_true * rho1
        r2 = R_atm_true + T_scat_true * rho2

        R_atm, T_scat, s_alb = extract_three_albedo(r0, r1, r2, rho1, rho2)

        self.assertAlmostEqual(R_atm, R_atm_true, places=10)
        self.assertAlmostEqual(T_scat, T_scat_true, places=6)
        self.assertAlmostEqual(s_alb, s_true, places=6)

    def test_array_input(self):
        """Extraction should work element-wise on arrays."""
        R_atm_true = np.array([0.02, 0.05, 0.10])
        T_scat_true = np.array([0.90, 0.80, 0.60])
        s_true = np.array([0.05, 0.10, 0.15])

        rho1, rho2 = 0.1, 0.5

        r0 = R_atm_true
        r1 = R_atm_true + T_scat_true * rho1 / (1 - s_true * rho1)
        r2 = R_atm_true + T_scat_true * rho2 / (1 - s_true * rho2)

        R_atm, T_scat, s_alb = extract_three_albedo(r0, r1, r2, rho1, rho2)

        np.testing.assert_allclose(R_atm, R_atm_true, atol=1e-10)
        np.testing.assert_allclose(T_scat, T_scat_true, atol=1e-6)
        np.testing.assert_allclose(s_alb, s_true, atol=1e-6)

    def test_high_path_radiance(self):
        """High R_atm (blue band, high AOD) should be recovered."""
        R_atm_true = 0.20
        T_scat_true = 0.50
        s_true = 0.20

        rho1, rho2 = 0.1, 0.5

        r0 = R_atm_true
        r1 = R_atm_true + T_scat_true * rho1 / (1 - s_true * rho1)
        r2 = R_atm_true + T_scat_true * rho2 / (1 - s_true * rho2)

        R_atm, T_scat, s_alb = extract_three_albedo(r0, r1, r2, rho1, rho2)

        self.assertAlmostEqual(R_atm, R_atm_true, places=10)
        self.assertAlmostEqual(T_scat, T_scat_true, places=8)
        self.assertAlmostEqual(s_alb, s_true, places=8)


def _make_synthetic_lut():
    """Create a synthetic LUT with known analytical values for testing.

    Uses simple analytical formulas:
        R_atm(wl, aod) = 0.01 * (550/wl)^4 + 0.1 * aod * (550/wl)^1.3
        T_scat(wl, aod) = exp(-0.5 * aod * (550/wl)^1.3)
        s(wl, aod) = 0.05 * (550/wl)^2 + 0.1 * aod
    """
    wavelengths = np.arange(400, 2501, 50, dtype=float)
    aod_grid = np.array([0.001, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5])
    n_wl = len(wavelengths)
    n_aod = len(aod_grid)

    R_atm = np.zeros((n_wl, n_aod))
    T_scat = np.zeros((n_wl, n_aod))
    s_alb = np.zeros((n_wl, n_aod))

    for i, wl in enumerate(wavelengths):
        for j, aod in enumerate(aod_grid):
            ratio = 550.0 / wl
            R_atm[i, j] = 0.01 * ratio**4 + 0.1 * aod * ratio**1.3
            T_scat[i, j] = np.exp(-0.5 * aod * ratio**1.3)
            s_alb[i, j] = min(1.0, 0.05 * ratio**2 + 0.1 * aod)

    return AtmosphericLUT(wavelengths, aod_grid, R_atm, T_scat, s_alb)


class TestInterpolationExact(unittest.TestCase):
    """At grid points, interpolation should return exact values."""

    def setUp(self):
        self.lut = _make_synthetic_lut()

    def test_exact_grid_point(self):
        """Interpolation at an exact grid point returns the stored value."""
        # Pick a grid point
        wl = 550.0
        aod = 0.2
        wl_idx = np.searchsorted(self.lut.wavelengths, wl)
        aod_idx = np.searchsorted(self.lut.aod_grid, aod)

        R, T, s = self.lut.interpolate(wl, aod)

        self.assertAlmostEqual(R, self.lut.R_atm[wl_idx, aod_idx], places=10)
        self.assertAlmostEqual(T, self.lut.T_scat[wl_idx, aod_idx], places=10)
        self.assertAlmostEqual(s, self.lut.s_alb[wl_idx, aod_idx], places=10)

    def test_multiple_grid_points(self):
        """Check several grid points."""
        for wl in [400.0, 800.0, 1200.0, 2000.0]:
            for aod in [0.001, 0.1, 0.5, 1.5]:
                R, T, s = self.lut.interpolate(wl, aod)
                wl_idx = np.searchsorted(self.lut.wavelengths, wl)
                aod_idx = np.searchsorted(self.lut.aod_grid, aod)
                self.assertAlmostEqual(
                    R, self.lut.R_atm[wl_idx, aod_idx], places=8,
                    msg=f"R_atm mismatch at wl={wl}, aod={aod}"
                )


class TestInterpolationBetween(unittest.TestCase):
    """Between grid points, results should be bounded by neighbors."""

    def setUp(self):
        self.lut = _make_synthetic_lut()

    def test_between_aod_values(self):
        """Interpolated R_atm at AOD=0.15 should be between values at 0.1 and 0.2."""
        wl = 550.0
        R_lo, _, _ = self.lut.interpolate(wl, 0.1)
        R_hi, _, _ = self.lut.interpolate(wl, 0.2)
        R_mid, _, _ = self.lut.interpolate(wl, 0.15)

        self.assertGreaterEqual(R_mid, min(R_lo, R_hi) - 1e-10)
        self.assertLessEqual(R_mid, max(R_lo, R_hi) + 1e-10)

    def test_between_wavelengths(self):
        """Interpolated at wl=575 should be between values at 550 and 600."""
        aod = 0.2
        R_lo, _, _ = self.lut.interpolate(550.0, aod)
        R_hi, _, _ = self.lut.interpolate(600.0, aod)
        R_mid, _, _ = self.lut.interpolate(575.0, aod)

        self.assertGreaterEqual(R_mid, min(R_lo, R_hi) - 1e-10)
        self.assertLessEqual(R_mid, max(R_lo, R_hi) + 1e-10)

    def test_array_aod_input(self):
        """Interpolation with a 2D AOD array should return same-shape output."""
        aod_arr = np.array([[0.1, 0.2], [0.3, 0.5]])
        R, T, s = self.lut.interpolate(550.0, aod_arr)
        self.assertEqual(R.shape, (2, 2))
        self.assertEqual(T.shape, (2, 2))
        self.assertEqual(s.shape, (2, 2))


class TestZeroAod(unittest.TestCase):
    """At AOD ~ 0, scattering is Rayleigh-only."""

    def setUp(self):
        self.lut = _make_synthetic_lut()

    def test_low_aod_small_ratm(self):
        """At near-zero AOD, R_atm should be small (Rayleigh only)."""
        R, T, s = self.lut.interpolate(850.0, 0.001)
        self.assertLess(R, 0.05)  # NIR Rayleigh is very small

    def test_low_aod_high_tscat(self):
        """At near-zero AOD, T_scat should be close to 1."""
        R, T, s = self.lut.interpolate(850.0, 0.001)
        self.assertGreater(T, 0.95)

    def test_low_aod_small_s(self):
        """At near-zero AOD, spherical albedo should be small."""
        R, T, s = self.lut.interpolate(850.0, 0.001)
        self.assertLess(s, 0.10)


class TestLutInversionIdentity(unittest.TestCase):
    """Forward(surface->TOA) then invert should recover original reflectance."""

    def setUp(self):
        self.lut = _make_synthetic_lut()

    def test_roundtrip_scalar(self):
        """Single-pixel roundtrip: forward then inverse."""
        wl = 650.0
        aod = 0.2
        r_surf_true = 0.15
        tg = 0.95  # Assumed gas transmission

        R_atm, T_scat, s = self.lut.interpolate(wl, aod)

        # Forward: r_toa = R_atm * tg + tg * T_scat * r_surf / (1 - s * r_surf)
        r_toa = R_atm * tg + tg * T_scat * r_surf_true / (1 - s * r_surf_true)

        # Inverse: r_surf = (r_toa - R_atm * tg) / (tg * T_scat + (r_toa - R_atm * tg) * s)
        numerator = r_toa - R_atm * tg
        denominator = tg * T_scat + numerator * s
        r_surf_recovered = numerator / denominator

        self.assertAlmostEqual(r_surf_recovered, r_surf_true, places=10)

    def test_roundtrip_array(self):
        """Array roundtrip at multiple AOD values."""
        wl = 550.0
        aod_arr = np.array([0.05, 0.1, 0.3, 0.5, 1.0])
        r_surf_true = 0.20
        tg = 0.90

        R_atm, T_scat, s = self.lut.interpolate(wl, aod_arr)

        r_toa = R_atm * tg + tg * T_scat * r_surf_true / (1 - s * r_surf_true)

        numerator = r_toa - R_atm * tg
        denominator = tg * T_scat + numerator * s
        r_surf_recovered = numerator / denominator

        np.testing.assert_allclose(r_surf_recovered, r_surf_true, atol=1e-10)

    def test_roundtrip_multiple_surfaces(self):
        """Roundtrip with different surface reflectance values."""
        wl = 480.0
        aod = 0.3
        surfaces = np.array([0.02, 0.05, 0.10, 0.20, 0.40])
        tg = 0.92

        R_atm, T_scat, s = self.lut.interpolate(wl, aod)

        r_toa = R_atm * tg + tg * T_scat * surfaces / (1 - s * surfaces)

        numerator = r_toa - R_atm * tg
        denominator = tg * T_scat + numerator * s
        r_surf_recovered = numerator / denominator

        np.testing.assert_allclose(r_surf_recovered, surfaces, atol=1e-10)


class TestSaveLoadRoundtrip(unittest.TestCase):
    """Save + load should preserve all LUT data exactly."""

    def test_roundtrip(self):
        lut = _make_synthetic_lut()

        with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
            tmp_path = f.name

        try:
            lut.save(tmp_path)
            lut2 = AtmosphericLUT.load(tmp_path)

            np.testing.assert_array_equal(lut.wavelengths, lut2.wavelengths)
            np.testing.assert_array_equal(lut.aod_grid, lut2.aod_grid)
            np.testing.assert_array_equal(lut.R_atm, lut2.R_atm)
            np.testing.assert_array_equal(lut.T_scat, lut2.T_scat)
            np.testing.assert_array_equal(lut.s_alb, lut2.s_alb)
        finally:
            os.unlink(tmp_path)

    def test_interpolation_after_load(self):
        """Interpolation results should be identical before and after save/load."""
        lut = _make_synthetic_lut()

        with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
            tmp_path = f.name

        try:
            lut.save(tmp_path)
            lut2 = AtmosphericLUT.load(tmp_path)

            for wl in [450.0, 700.0, 1500.0]:
                for aod in [0.05, 0.3, 1.0]:
                    R1, T1, s1 = lut.interpolate(wl, aod)
                    R2, T2, s2 = lut2.interpolate(wl, aod)
                    self.assertAlmostEqual(R1, R2, places=14)
                    self.assertAlmostEqual(T1, T2, places=14)
                    self.assertAlmostEqual(s1, s2, places=14)
        finally:
            os.unlink(tmp_path)


class TestCacheKey(unittest.TestCase):
    """Cache key generation should be deterministic and unique."""

    def test_deterministic(self):
        key1 = AtmosphericLUT.cache_key(30, 0, 45, 1013, 'continental',
                                         350, 2500, 10)
        key2 = AtmosphericLUT.cache_key(30, 0, 45, 1013, 'continental',
                                         350, 2500, 10)
        self.assertEqual(key1, key2)

    def test_different_params_different_key(self):
        key1 = AtmosphericLUT.cache_key(30, 0, 45, 1013, 'continental',
                                         350, 2500, 10)
        key2 = AtmosphericLUT.cache_key(40, 0, 45, 1013, 'continental',
                                         350, 2500, 10)
        self.assertNotEqual(key1, key2)


class TestPhysicalBounds(unittest.TestCase):
    """Interpolated values should always be within physical bounds."""

    def setUp(self):
        self.lut = _make_synthetic_lut()

    def test_bounds_sweep(self):
        """Sweep many (wl, aod) pairs and check physical bounds."""
        for wl in np.arange(400, 2501, 100):
            for aod in [0.001, 0.05, 0.1, 0.3, 0.5, 1.0, 1.5]:
                R, T, s = self.lut.interpolate(float(wl), aod)
                self.assertGreaterEqual(R, 0.0,
                    msg=f"R_atm < 0 at wl={wl}, aod={aod}")
                self.assertGreaterEqual(T, 0.0,
                    msg=f"T_scat < 0 at wl={wl}, aod={aod}")
                self.assertLessEqual(T, 1.0,
                    msg=f"T_scat > 1 at wl={wl}, aod={aod}")
                self.assertGreaterEqual(s, 0.0,
                    msg=f"s_alb < 0 at wl={wl}, aod={aod}")
                self.assertLessEqual(s, 1.0,
                    msg=f"s_alb > 1 at wl={wl}, aod={aod}")

    def test_clamped_outside_range(self):
        """Values outside the grid should be clamped, not extrapolated."""
        # Below minimum wavelength
        R, T, s = self.lut.interpolate(300.0, 0.1)
        self.assertGreaterEqual(R, 0.0)
        self.assertGreaterEqual(T, 0.0)
        self.assertLessEqual(T, 1.0)

        # Above maximum AOD
        R, T, s = self.lut.interpolate(550.0, 3.0)
        self.assertGreaterEqual(R, 0.0)
        self.assertLessEqual(s, 1.0)


if __name__ == '__main__':
    unittest.main()
