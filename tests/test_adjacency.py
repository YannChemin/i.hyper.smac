#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for adjacency effect correction (Improvement 5).
"""

import unittest
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from adjacency import (
    compute_environmental_reflectance,
    compute_direct_transmittance,
    adjacency_correction,
)


class TestComputeEnvironmentalReflectance(unittest.TestCase):
    """Tests for compute_environmental_reflectance()."""

    def test_uniform_scene(self):
        """Uniform reflectance should give the same environmental reflectance."""
        r = np.full((50, 50), 0.30)
        r_env = compute_environmental_reflectance(r, pixel_size=30.0,
                                                   psf_radius_km=1.0)
        # Interior pixels should be very close to 0.30
        np.testing.assert_allclose(r_env[10:40, 10:40], 0.30, atol=0.01)

    def test_bright_patch(self):
        """Environmental reflectance of a bright patch should be lower than
        the patch itself (diluted by dark surroundings)."""
        r = np.full((100, 100), 0.05)
        r[40:60, 40:60] = 0.50  # bright patch
        r_env = compute_environmental_reflectance(r, pixel_size=30.0,
                                                   psf_radius_km=1.0)
        # Center of bright patch: env should be < 0.50 due to dark neighbours
        self.assertLess(r_env[50, 50], 0.50)
        # But still higher than pure dark background
        self.assertGreater(r_env[50, 50], 0.05)


class TestComputeDirectTransmittance(unittest.TestCase):
    """Tests for compute_direct_transmittance()."""

    def test_clear_atmosphere_high_transmittance(self):
        """Low AOD at NIR should give near-unity direct transmittance."""
        T = compute_direct_transmittance(850, aod550=0.05,
                                          pressure=1013.25,
                                          sza=30, vza=0)
        self.assertGreater(T, 0.85)

    def test_T_dir_increases_with_wavelength(self):
        """Less scattering at longer wavelengths → higher T_dir."""
        T_400 = compute_direct_transmittance(400, aod550=0.2,
                                              pressure=1013.25,
                                              sza=30, vza=0)
        T_800 = compute_direct_transmittance(800, aod550=0.2,
                                              pressure=1013.25,
                                              sza=30, vza=0)
        self.assertGreater(float(T_800), float(T_400))

    def test_zero_aod_rayleigh_only(self):
        """With zero AOD, transmittance should be Rayleigh-only."""
        T = compute_direct_transmittance(550, aod550=0.0,
                                          pressure=1013.25,
                                          sza=30, vza=0)
        # Should be high but not exactly 1 due to Rayleigh
        self.assertGreater(T, 0.5)
        self.assertLess(T, 1.0)


class TestAdjacencyCorrection(unittest.TestCase):
    """Tests for adjacency_correction()."""

    def test_uniform_scene_no_change(self):
        """Uniform reflectance should produce negligible correction."""
        r = np.full((50, 50), 0.20)
        T_scat = 0.80
        s = 0.05

        r_adj = adjacency_correction(r, T_scat, s, wavelength=550,
                                      aod550=0.2, pressure=1013.25,
                                      sza=30, vza=0,
                                      pixel_size=30.0, psf_radius_km=1.0)
        # Should be nearly unchanged
        np.testing.assert_allclose(r_adj[10:40, 10:40], 0.20, atol=0.01)

    def test_bright_patch_corrected(self):
        """Bright patch in dark background should have correction applied."""
        r = np.full((100, 100), 0.05)
        r[40:60, 40:60] = 0.40
        T_scat = 0.70
        s = 0.08

        r_adj = adjacency_correction(r, T_scat, s, wavelength=550,
                                      aod550=0.3, pressure=1013.25,
                                      sza=30, vza=0,
                                      pixel_size=30.0, psf_radius_km=1.0)

        # Center of bright patch: r_boa > r_env, so correction is positive
        # → corrected reflectance should be higher than original
        self.assertGreater(r_adj[50, 50], r[50, 50])

    def test_clear_atmosphere_minimal(self):
        """Near-zero AOD at NIR should produce T_dir ≈ T_scat, minimal T_diff."""
        r = np.full((50, 50), 0.10)
        r[20:30, 20:30] = 0.50
        T_scat = 0.95
        s = 0.02

        r_adj = adjacency_correction(r, T_scat, s, wavelength=850,
                                      aod550=0.01, pressure=1013.25,
                                      sza=30, vza=0,
                                      pixel_size=30.0, psf_radius_km=1.0)

        # Correction should be very small
        max_diff = np.nanmax(np.abs(r_adj - r))
        self.assertLess(max_diff, 0.05)

    def test_nan_preserved(self):
        """NaN values in the input should be preserved."""
        r = np.full((20, 20), 0.20)
        r[5, 5] = np.nan
        T_scat = 0.80
        s = 0.05

        r_adj = adjacency_correction(r, T_scat, s, wavelength=550,
                                      aod550=0.2, pressure=1013.25,
                                      sza=30, vza=0,
                                      pixel_size=30.0, psf_radius_km=1.0)
        self.assertTrue(np.isnan(r_adj[5, 5]))


if __name__ == '__main__':
    unittest.main()
