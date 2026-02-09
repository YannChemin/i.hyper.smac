#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for atmospheric correction helper functions.

Tests the transmission masking, reflectance clamping, AOD validation,
and blue-band AOD tapering introduced to fix systematic SMAC issues.
"""

import unittest
import sys
import os
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from smac import smac_inv


# ---------------------------------------------------------------------------
# Re-implement the helpers here so tests don't depend on GRASS module import.
# These must stay in sync with the definitions in i.hyper.smac.py.
# ---------------------------------------------------------------------------

def compute_band_transmission(coefs, sza, vza, pressure, aod550, wvc, o3):
    us = np.cos(np.radians(sza))
    uv = np.cos(np.radians(vza))
    Peq = pressure / 1013.25
    m = 1.0 / us + 1.0 / uv

    th2o = np.exp(coefs.ah2o * ((wvc * m) ** coefs.nh2o))
    to3 = np.exp(coefs.ao3 * ((o3 * m) ** coefs.no3))

    uo2 = Peq ** coefs.po2
    to2 = np.exp(coefs.ao2 * ((uo2 * m) ** coefs.no2))
    uco2 = Peq ** coefs.pco2
    tco2 = np.exp(coefs.aco2 * ((uco2 * m) ** coefs.nco2))
    uch4 = Peq ** coefs.pch4
    tch4 = np.exp(coefs.ach4 * ((uch4 * m) ** coefs.nch4))
    uno2 = Peq ** coefs.pno2
    tno2 = np.exp(coefs.ano2 * ((uno2 * m) ** coefs.nno2))
    uco = Peq ** coefs.pco
    tco = np.exp(coefs.aco * ((uco * m) ** coefs.nco))

    tg = th2o * to3 * to2 * tco2 * tch4 * tco * tno2

    ttetas = (coefs.a0T + coefs.a1T * aod550 / us
              + (coefs.a2T * Peq + coefs.a3T) / (1.0 + us))
    ttetav = (coefs.a0T + coefs.a1T * aod550 / uv
              + (coefs.a2T * Peq + coefs.a3T) / (1.0 + uv))

    return float(tg * ttetas * ttetav)


def compute_blue_aod_taper(wavelength, aod):
    if wavelength >= 500.0:
        return aod
    fraction = max(0.7, 0.7 + 0.3 * (wavelength - 400.0) / 100.0)
    return aod * fraction


AOD_MAX = 1.5
TRANSMISSION_THRESHOLD = 0.10


# ---------------------------------------------------------------------------
# Mock coefficient objects
# ---------------------------------------------------------------------------

class _CoefBase:
    """Transparent-window SMAC coefficients (no gaseous absorption)."""
    def __init__(self):
        self.ah2o = 0.0;  self.nh2o = 0.5
        self.ao3 = 0.0;   self.no3 = 1.0
        self.ao2 = 0.0;   self.no2 = 0.0;  self.po2 = 1.0
        self.aco2 = 0.0;  self.nco2 = 0.0; self.pco2 = 1.0
        self.ach4 = 0.0;  self.nch4 = 0.0; self.pch4 = 1.0
        self.ano2 = 0.0;  self.nno2 = 0.0; self.pno2 = 1.0
        self.aco = 0.0;   self.nco = 0.0;  self.pco = 1.0
        self.a0s = 0.05;  self.a1s = 0.01; self.a2s = 0.001; self.a3s = 0.0
        self.a0T = 1.0;   self.a1T = -0.15; self.a2T = 0.0; self.a3T = -0.05
        self.taur = 0.05; self.sr = 0.05
        self.a0taup = 0.0; self.a1taup = 1.0
        self.wo = 0.9;    self.gc = 0.7
        self.a0P = 0.1;   self.a1P = 0.0; self.a2P = 0.0
        self.a3P = 0.0;   self.a4P = 0.0
        self.Rest1 = 0.0; self.Rest2 = 0.0; self.Rest3 = 0.0; self.Rest4 = 0.0
        self.Resr1 = 0.0; self.Resr2 = 0.0; self.Resr3 = 0.0
        self.Resa1 = 0.0; self.Resa2 = 0.0; self.Resa3 = 0.0; self.Resa4 = 0.0


class _CoefOpaque(_CoefBase):
    """Coefficients mimicking a strong H2O absorption band."""
    def __init__(self):
        super().__init__()
        # Strong H2O absorption → very low gas transmission
        self.ah2o = -2.5
        self.nh2o = 0.6


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestComputeBandTransmission(unittest.TestCase):
    """Tests for compute_band_transmission()."""

    def test_transparent_window(self):
        """A clean atmospheric window should have transmission near 1."""
        coefs = _CoefBase()
        T = compute_band_transmission(coefs, sza=30, vza=0,
                                      pressure=1013.25, aod550=0.1,
                                      wvc=2.0, o3=0.3)
        self.assertGreater(T, 0.5)
        self.assertLessEqual(T, 1.0)

    def test_opaque_band(self):
        """A deep absorption band should have very low transmission."""
        coefs = _CoefOpaque()
        T = compute_band_transmission(coefs, sza=40, vza=30,
                                      pressure=1013.25, aod550=0.15,
                                      wvc=3.0, o3=0.3)
        self.assertLess(T, TRANSMISSION_THRESHOLD)

    def test_high_aod_lowers_transmission(self):
        """Increasing AOD should decrease scattering transmission."""
        coefs = _CoefBase()
        T_low = compute_band_transmission(coefs, sza=30, vza=0,
                                          pressure=1013.25, aod550=0.05,
                                          wvc=2.0, o3=0.3)
        T_high = compute_band_transmission(coefs, sza=30, vza=0,
                                           pressure=1013.25, aod550=0.8,
                                           wvc=2.0, o3=0.3)
        self.assertGreater(T_low, T_high)

    def test_returns_scalar(self):
        """Result should be a Python float, not an array."""
        coefs = _CoefBase()
        T = compute_band_transmission(coefs, sza=30, vza=0,
                                      pressure=1013.25, aod550=0.1,
                                      wvc=2.0, o3=0.3)
        self.assertIsInstance(T, float)


class TestComputeBlueAodTaper(unittest.TestCase):
    """Tests for compute_blue_aod_taper()."""

    def test_no_taper_above_500(self):
        """No taper at 500 nm or above."""
        aod = 0.5
        self.assertEqual(compute_blue_aod_taper(500.0, aod), aod)
        self.assertEqual(compute_blue_aod_taper(600.0, aod), aod)
        self.assertEqual(compute_blue_aod_taper(1000.0, aod), aod)

    def test_taper_at_450(self):
        """At 450 nm, AOD should be reduced to 85%."""
        aod = 1.0
        result = compute_blue_aod_taper(450.0, aod)
        self.assertAlmostEqual(result, 0.85, places=5)

    def test_taper_at_400(self):
        """At 400 nm, AOD should be reduced to 70%."""
        aod = 1.0
        result = compute_blue_aod_taper(400.0, aod)
        self.assertAlmostEqual(result, 0.70, places=5)

    def test_taper_capped_below_400(self):
        """Below 400 nm, AOD fraction should stay at 70%."""
        aod = 1.0
        result = compute_blue_aod_taper(350.0, aod)
        self.assertAlmostEqual(result, 0.70, places=5)

    def test_zero_aod(self):
        """Zero AOD should stay zero regardless of wavelength."""
        self.assertEqual(compute_blue_aod_taper(420.0, 0.0), 0.0)


class TestAodClamping(unittest.TestCase):
    """Tests for AOD upper-bound validation (Fix 3)."""

    def test_normal_aod_unchanged(self):
        """AOD within bounds should not be modified."""
        aod = 0.3
        self.assertEqual(min(aod, AOD_MAX), aod)

    def test_extreme_aod_clamped(self):
        """AOD above maximum should be clamped."""
        aod = 3.5
        self.assertEqual(min(aod, AOD_MAX), AOD_MAX)

    def test_boundary_aod(self):
        """AOD exactly at the limit should pass through."""
        self.assertEqual(min(AOD_MAX, AOD_MAX), AOD_MAX)


class TestReflectanceClamping(unittest.TestCase):
    """Tests for post-SMAC reflectance clamping (Fix 2)."""

    def test_normal_values_unchanged(self):
        """Values within [-0.01, 1.5] should not be clipped."""
        data = np.array([0.0, 0.1, 0.5, 1.0])
        clipped = np.clip(data, -0.01, 1.5)
        np.testing.assert_array_equal(data, clipped)

    def test_negative_clipped(self):
        """Strongly negative reflectance should be clipped to -0.01."""
        data = np.array([-5.0, -0.02, -11.0])
        clipped = np.clip(data, -0.01, 1.5)
        self.assertTrue(np.all(clipped >= -0.01))

    def test_high_values_clipped(self):
        """Reflectance above 1.5 should be clipped."""
        data = np.array([2.0, 3.5, 10.0])
        clipped = np.clip(data, -0.01, 1.5)
        self.assertTrue(np.all(clipped <= 1.5))

    def test_small_negative_preserved(self):
        """Small negative values (adjacency/noise) should be preserved."""
        data = np.array([-0.005, -0.009])
        clipped = np.clip(data, -0.01, 1.5)
        np.testing.assert_array_almost_equal(data, clipped)


class TestTransmissionMaskingWithSMAC(unittest.TestCase):
    """Integration test: verify that low-transmission bands produce NaN."""

    def test_opaque_band_gives_nan(self):
        """Simulates the band-loop logic: opaque band → NaN output."""
        coefs = _CoefOpaque()
        sza, vza = 40.0, 30.0
        pressure, aod, wvc, o3 = 1013.25, 0.15, 3.0, 0.3

        T = compute_band_transmission(coefs, sza, vza, pressure, aod, wvc, o3)

        r_toa = np.array([0.1, 0.2, 0.3])
        if T < TRANSMISSION_THRESHOLD:
            refl_boa = np.full_like(r_toa, np.nan)
        else:
            refl_boa = smac_inv(r_toa, sza, 180, vza, 0,
                                pressure, aod, o3, wvc, coefs)
            refl_boa = np.clip(refl_boa, -0.01, 1.5)

        self.assertTrue(np.all(np.isnan(refl_boa)))

    def test_transparent_band_gives_valid(self):
        """Simulates the band-loop logic: transparent band → valid output."""
        coefs = _CoefBase()
        sza, vza = 30.0, 0.0
        pressure, aod, wvc, o3 = 1013.25, 0.1, 2.0, 0.3

        T = compute_band_transmission(coefs, sza, vza, pressure, aod, wvc, o3)

        r_toa = np.array([0.15, 0.20, 0.25])
        if T < TRANSMISSION_THRESHOLD:
            refl_boa = np.full_like(r_toa, np.nan)
        else:
            aod_band = compute_blue_aod_taper(650.0, aod)
            refl_boa = smac_inv(r_toa, sza, 180, vza, 0,
                                pressure, aod_band, o3, wvc, coefs)
            refl_boa = np.clip(refl_boa, -0.01, 1.5)

        self.assertFalse(np.any(np.isnan(refl_boa)))
        self.assertTrue(np.all(refl_boa >= -0.01))
        self.assertTrue(np.all(refl_boa <= 1.5))


if __name__ == '__main__':
    unittest.main()
