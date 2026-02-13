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

    return tg * ttetas * ttetav


def compute_blue_aod_taper(wavelength, aod):
    if wavelength >= 650.0:
        return aod
    alpha = np.minimum(2.0, 2.0 * (650.0 - wavelength) / 250.0)
    return aod / (1.0 + alpha * aod)


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
        """Result should be a scalar (0-d) when all inputs are scalar."""
        coefs = _CoefBase()
        T = compute_band_transmission(coefs, sza=30, vza=0,
                                      pressure=1013.25, aod550=0.1,
                                      wvc=2.0, o3=0.3)
        self.assertEqual(np.ndim(T), 0)


class TestComputeBlueAodTaper(unittest.TestCase):
    """Tests for compute_blue_aod_taper() saturation model."""

    def test_no_taper_above_650(self):
        """No taper at 650 nm or above."""
        aod = 0.5
        self.assertEqual(compute_blue_aod_taper(650.0, aod), aod)
        self.assertEqual(compute_blue_aod_taper(700.0, aod), aod)
        self.assertEqual(compute_blue_aod_taper(1000.0, aod), aod)

    def test_low_aod_minimal_change(self):
        """At low AOD (0.1), correction should be small (<20% reduction)."""
        aod = 0.1
        result = compute_blue_aod_taper(400.0, aod)
        # alpha=2.0 at 400nm → 0.1/(1+0.2) = 0.0833
        self.assertAlmostEqual(result, 0.0833, places=3)
        self.assertGreater(result, aod * 0.8)  # less than 20% reduction

    def test_high_aod_strong_reduction(self):
        """At high AOD (0.664), correction should exceed 50% at 400nm."""
        aod = 0.664
        result = compute_blue_aod_taper(400.0, aod)
        # alpha=2.0 → 0.664/(1+1.328) = 0.285
        self.assertAlmostEqual(result, 0.285, places=2)
        self.assertLess(result, aod * 0.5)  # more than 50% reduction

    def test_saturation_at_500nm(self):
        """At 500nm with AOD=1.0, alpha=1.2 → 1/(1+1.2) = 0.4545."""
        result = compute_blue_aod_taper(500.0, 1.0)
        self.assertAlmostEqual(result, 1.0 / 2.2, places=4)

    def test_zero_aod(self):
        """Zero AOD should stay zero regardless of wavelength."""
        self.assertEqual(compute_blue_aod_taper(420.0, 0.0), 0.0)

    def test_monotonic_in_wavelength(self):
        """For fixed AOD, effective AOD should increase with wavelength."""
        aod = 0.5
        wavelengths = [400, 450, 500, 550, 600, 650]
        results = [compute_blue_aod_taper(wl, aod) for wl in wavelengths]
        for i in range(len(results) - 1):
            self.assertLessEqual(results[i], results[i + 1])


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


class TestSmacInvArray(unittest.TestCase):
    """Tests for smac_inv() with 2D array atmospheric parameters."""

    def test_array_pressure_and_aod(self):
        """smac_inv should handle 2D pressure and AOD arrays."""
        coefs = _CoefBase()
        r_toa = np.full((3, 4), 0.2)
        pressure = np.full((3, 4), 1013.25)
        aod = np.full((3, 4), 0.15)
        wvc = np.full((3, 4), 2.0)
        r_surf = smac_inv(r_toa, 30, 180, 0, 0, pressure, aod, 0.3, wvc, coefs)
        self.assertEqual(r_surf.shape, (3, 4))
        self.assertTrue(np.all(np.isfinite(r_surf)))

    def test_uniform_array_matches_scalar(self):
        """Uniform-array inputs should produce the same result as scalars."""
        coefs = _CoefBase()
        r_toa = np.array([0.15, 0.20, 0.25])
        sza, vza = 30.0, 0.0
        pres, aod_val, o3_val, wvc_val = 1013.25, 0.1, 0.3, 2.0

        r_scalar = smac_inv(r_toa, sza, 180, vza, 0,
                            pres, aod_val, o3_val, wvc_val, coefs)

        shape = (3,)
        r_array = smac_inv(r_toa, sza, 180, vza, 0,
                           np.full(shape, pres),
                           np.full(shape, aod_val),
                           o3_val,
                           np.full(shape, wvc_val),
                           coefs)
        np.testing.assert_allclose(r_array, r_scalar, rtol=1e-12)


class TestComputeBandTransmissionArray(unittest.TestCase):
    """Tests for compute_band_transmission() with 2D array inputs."""

    def test_array_pressure_returns_2d(self):
        """Passing a 2D pressure array should return a 2D result."""
        coefs = _CoefBase()
        pressure = np.array([[1013.25, 900.0], [850.0, 1013.25]])
        T = compute_band_transmission(coefs, sza=30, vza=0,
                                      pressure=pressure, aod550=0.1,
                                      wvc=2.0, o3=0.3)
        self.assertEqual(T.shape, (2, 2))

    def test_array_aod_returns_2d(self):
        """Passing a 2D AOD array should return a 2D result."""
        coefs = _CoefBase()
        aod = np.array([[0.05, 0.1], [0.3, 0.5]])
        T = compute_band_transmission(coefs, sza=30, vza=0,
                                      pressure=1013.25, aod550=aod,
                                      wvc=2.0, o3=0.3)
        self.assertEqual(T.shape, (2, 2))

    def test_uniform_array_matches_scalar(self):
        """A uniform 2D array should give the same result as the scalar."""
        coefs = _CoefBase()
        sza, vza = 30.0, 10.0
        pres, aod_val, wvc_val, o3_val = 1013.25, 0.15, 2.5, 0.3

        T_scalar = compute_band_transmission(coefs, sza, vza,
                                             pres, aod_val, wvc_val, o3_val)

        shape = (3, 4)
        T_array = compute_band_transmission(
            coefs, sza, vza,
            np.full(shape, pres),
            np.full(shape, aod_val),
            np.full(shape, wvc_val),
            np.full(shape, o3_val),
        )
        self.assertEqual(T_array.shape, shape)
        np.testing.assert_allclose(T_array, float(T_scalar), rtol=1e-12)


class TestComputeBlueAodTaperArray(unittest.TestCase):
    """Tests for compute_blue_aod_taper() with array AOD."""

    def test_array_aod_returns_array(self):
        """Passing an array AOD should return an array of same shape."""
        aod = np.array([0.1, 0.3, 0.5])
        result = compute_blue_aod_taper(450.0, aod)
        self.assertEqual(result.shape, (3,))

    def test_no_taper_above_650_array(self):
        """Above 650 nm, array AOD should pass through unchanged."""
        aod = np.array([0.1, 0.5, 1.0])
        result = compute_blue_aod_taper(700.0, aod)
        np.testing.assert_array_equal(result, aod)

    def test_scalar_and_uniform_array_match(self):
        """Scalar and uniform-array AOD should give identical results."""
        aod_scalar = 0.3
        aod_array = np.full((2, 3), aod_scalar)
        wl = 420.0

        r_scalar = compute_blue_aod_taper(wl, aod_scalar)
        r_array = compute_blue_aod_taper(wl, aod_array)

        np.testing.assert_allclose(r_array, float(r_scalar), rtol=1e-12)


# ---------------------------------------------------------------------------
# Coupling correction helper (mirrors i.hyper.smac.py)
# ---------------------------------------------------------------------------

def compute_coupling_correction(wavelength, tg, aod550, pressure,
                                 aerosol_model='continental', k=0.07):
    wl_um = wavelength / 1000.0
    tau_r = 0.008569 * wl_um**(-4) * (1 + 0.0113 * wl_um**(-2))
    tau_r = tau_r * np.asarray(pressure) / 1013.25
    tau_a = np.asarray(aod550) * (wavelength / 550.0)**(-1.3)
    tau_scat = tau_r + tau_a
    return np.asarray(tg) ** (1.0 + k * tau_scat)


# ---------------------------------------------------------------------------
# Tests for coupling correction (Improvement 4)
# ---------------------------------------------------------------------------

class TestCouplingCorrection(unittest.TestCase):
    """Tests for compute_coupling_correction()."""

    def test_no_correction_at_zero_scattering(self):
        """tg_eff should equal tg when AOD=0 at long wavelengths (low Rayleigh)."""
        tg = 0.95
        # At 2000nm, Rayleigh is negligible; with aod=0, tau_scat ≈ 0
        tg_eff = compute_coupling_correction(2000.0, tg, aod550=0.0,
                                              pressure=1013.25)
        self.assertAlmostEqual(float(tg_eff), tg, places=4)

    def test_correction_increases_absorption(self):
        """tg_eff should be less than tg when scattering is present."""
        tg = 0.90
        tg_eff = compute_coupling_correction(500.0, tg, aod550=0.3,
                                              pressure=1013.25)
        self.assertLess(float(tg_eff), tg)

    def test_stronger_at_short_wavelengths(self):
        """Correction should be larger at 400nm than at 800nm (more Rayleigh)."""
        tg = 0.90
        tg_400 = compute_coupling_correction(400.0, tg, aod550=0.2,
                                              pressure=1013.25)
        tg_800 = compute_coupling_correction(800.0, tg, aod550=0.2,
                                              pressure=1013.25)
        # Stronger correction = lower tg_eff
        self.assertLess(float(tg_400), float(tg_800))

    def test_weak_absorption_minimal(self):
        """When tg is very close to 1, correction should be negligible."""
        tg = 0.999
        tg_eff = compute_coupling_correction(500.0, tg, aod550=0.3,
                                              pressure=1013.25)
        self.assertAlmostEqual(float(tg_eff), tg, places=3)

    def test_array_input(self):
        """Should work with 2D AOD/pressure arrays."""
        tg = 0.90
        aod = np.array([[0.1, 0.2], [0.3, 0.5]])
        pressure = np.array([[1013.25, 900.0], [850.0, 1013.25]])
        tg_eff = compute_coupling_correction(500.0, tg, aod, pressure)
        self.assertEqual(tg_eff.shape, (2, 2))
        self.assertTrue(np.all(tg_eff <= tg))
        self.assertTrue(np.all(tg_eff > 0))


# ---------------------------------------------------------------------------
# Tests for blue LUT hybrid (Improvement 2)
# ---------------------------------------------------------------------------

class TestBlueHybridConcept(unittest.TestCase):
    """Verify the hybrid LUT+SMAC gas approach produces reasonable results."""

    def test_lut_inversion_positive_for_dark_surface(self):
        """LUT-style inversion should give non-negative reflectance for
        typical dark vegetation TOA signal at blue wavelengths."""
        # Simulate: r_toa = R_atm * tg + tg * T_scat * rho / (1 - rho * s)
        # For dark veg at 450nm: rho_surf ≈ 0.03, heavy atmosphere
        R_atm = 0.15     # Path radiance (DISORT-like, reasonable for blue)
        T_scat = 0.60     # Two-way scattering transmittance
        s = 0.10           # Spherical albedo
        tg = 0.95          # Gas transmission
        rho_true = 0.03    # True surface reflectance

        # Simulate TOA reflectance
        r_toa = R_atm * tg + tg * T_scat * rho_true / (1.0 - rho_true * s)

        # Inversion
        numerator = r_toa - R_atm * tg
        denominator = tg * T_scat + numerator * s
        rho_recovered = numerator / denominator

        self.assertGreater(rho_recovered, -0.01)
        self.assertAlmostEqual(rho_recovered, rho_true, places=3)

    def test_blend_region_smooth(self):
        """A simple linear blend between LUT and SMAC should produce
        intermediate values between 550-650nm."""
        val_lut = 0.05      # Blue side (LUT)
        val_smac = 0.06     # Red side (SMAC)
        blend_start, blend_end = 550, 650

        for wl in [550, 575, 600, 625, 650]:
            if wl <= blend_start:
                blended = val_lut
            elif wl >= blend_end:
                blended = val_smac
            else:
                f = (wl - blend_start) / (blend_end - blend_start)
                blended = (1 - f) * val_lut + f * val_smac

            self.assertGreaterEqual(blended, min(val_lut, val_smac))
            self.assertLessEqual(blended, max(val_lut, val_smac))


if __name__ == '__main__':
    unittest.main()
