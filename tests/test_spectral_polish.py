#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for spectral polishing (Improvement 3).
"""

import unittest
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from spectral_polish import (
    compute_quality_weights,
    spectral_polish,
    interpolate_flagged,
)


class TestComputeQualityWeights(unittest.TestCase):
    """Tests for compute_quality_weights()."""

    def test_full_transmission_high_quality(self):
        """Bands with tg=1 should have quality near 1."""
        wl = np.arange(400, 900, 10.0)
        tg = np.ones_like(wl)
        qw = compute_quality_weights(wl, tg)
        self.assertTrue(np.all(qw > 0.9))

    def test_zero_transmission_zero_quality(self):
        """Bands with tg=0 should have quality 0."""
        wl = np.arange(400, 900, 10.0)
        tg = np.zeros_like(wl)
        qw = compute_quality_weights(wl, tg)
        self.assertTrue(np.all(qw == 0.0))

    def test_edge_penalty_reduces_quality(self):
        """Rapid change in tg (absorption edge) should reduce quality."""
        wl = np.arange(750, 780, 1.0)
        tg = np.ones_like(wl)
        # Simulate O2-A absorption edge: sharp drop at 762nm
        tg[wl >= 762] = 0.3

        qw = compute_quality_weights(wl, tg)
        # The band at the edge (762nm) should have lower quality than 750nm
        idx_750 = 0
        idx_edge = np.argmin(np.abs(wl - 762))
        self.assertLess(qw[idx_edge], qw[idx_750])


class TestSpectralPolish(unittest.TestCase):
    """Tests for spectral_polish()."""

    def test_clean_spectrum_unchanged(self):
        """A smooth vegetation-like spectrum should not be modified."""
        wl = np.arange(400, 900, 10.0)
        refl = 0.05 + 0.30 * (1 / (1 + np.exp(-(wl - 700) / 20)))
        polished, flags = spectral_polish(refl, wl)
        self.assertFalse(np.any(flags))
        np.testing.assert_array_almost_equal(polished, refl, decimal=6)

    def test_single_spike_detected(self):
        """A single spike (0.73 among 0.35 neighbours) should be flagged."""
        wl = np.arange(750, 780, 1.0)
        refl = np.full_like(wl, 0.35)
        refl[16] = 0.73  # spike at 766nm

        _, flags = spectral_polish(refl, wl, mad_threshold=3.0, replace=False)
        self.assertTrue(flags[16])

    def test_spike_replaced_interpolated(self):
        """Flagged spike should be replaced with interpolated value."""
        wl = np.arange(750, 780, 1.0)
        refl = np.full_like(wl, 0.35)
        refl[16] = 0.73  # spike at 766nm

        polished, flags = spectral_polish(refl, wl, replace=True)
        if flags[16]:
            self.assertAlmostEqual(polished[16], 0.35, places=1)

    def test_nan_preserved(self):
        """NaN (opaque) bands should stay NaN after polishing."""
        wl = np.arange(400, 500, 5.0)
        refl = np.full_like(wl, 0.05)
        refl[5] = np.nan
        refl[6] = np.nan

        polished, flags = spectral_polish(refl, wl)
        self.assertTrue(np.isnan(polished[5]))
        self.assertTrue(np.isnan(polished[6]))

    def test_multiple_consecutive_spikes(self):
        """Multiple consecutive spikes should be handled gracefully."""
        wl = np.arange(400, 500, 2.0)
        refl = np.full_like(wl, 0.10)
        refl[10:13] = 0.50  # Three consecutive spikes

        polished, flags = spectral_polish(refl, wl, replace=True)
        # At least some of the spikes should be flagged
        self.assertTrue(np.any(flags[10:13]))
        # Polished values should be closer to 0.10 than original 0.50
        self.assertTrue(np.all(polished[10:13] < 0.40))

    def test_batch_2d_input(self):
        """[n_pixels, n_bands] shape should be preserved."""
        wl = np.arange(400, 500, 5.0)
        n_pixels = 10
        refl = np.full((n_pixels, len(wl)), 0.20)
        # Add a spike to each pixel at different positions
        for i in range(n_pixels):
            refl[i, i % len(wl)] = 0.80

        polished, flags = spectral_polish(refl, wl, replace=True)
        self.assertEqual(polished.shape, (n_pixels, len(wl)))
        self.assertEqual(flags.shape, (n_pixels, len(wl)))


class TestInterpolateFlagged(unittest.TestCase):
    """Tests for interpolate_flagged()."""

    def test_no_flags_no_change(self):
        """No flagged bands → no change."""
        wl = np.arange(400, 500, 10.0)
        refl = np.array([[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])
        flags = np.zeros_like(refl, dtype=bool)
        result = interpolate_flagged(refl, flags, wl)
        np.testing.assert_array_equal(result, refl)

    def test_single_flagged_interpolated(self):
        """A single flagged band should be linearly interpolated."""
        wl = np.array([400, 410, 420, 430, 440])
        refl = np.array([[0.1, 0.2, 999.0, 0.4, 0.5]])
        flags = np.array([[False, False, True, False, False]])
        result = interpolate_flagged(refl, flags, wl)
        # Should interpolate to ~0.3
        self.assertAlmostEqual(result[0, 2], 0.3, places=2)


if __name__ == '__main__':
    unittest.main()
