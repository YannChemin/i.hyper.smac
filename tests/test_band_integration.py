#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for band-integrated SMAC coefficients (Improvement 1).

Tests the band integration helper, FWHM-specific file lookup, and
comment-line handling in the coefficient parser.
"""

import unittest
import sys
import os
import tempfile
import shutil
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from smac_coef_generator import _band_integrate_result, SMACCoefficients


class TestBandIntegrateResult(unittest.TestCase):
    """Tests for _band_integrate_result() Gaussian SRF weighting."""

    def test_smooth_region_matches_center(self):
        """At 850nm over a flat spectrum, band-integrated ≈ center value."""
        wls = np.arange(840, 861, 1.0)
        values = np.full_like(wls, 100.0)  # Flat spectrum
        result = _band_integrate_result(850, 6.0, wls, values)
        self.assertAlmostEqual(result, 100.0, places=3)

    def test_linear_spectrum_center_weighted(self):
        """A linear ramp should be weighted toward the center wavelength."""
        wls = np.arange(840, 861, 1.0)
        values = wls * 1.0  # Linearly increasing
        result = _band_integrate_result(850, 6.0, wls, values)
        # Should be close to center value of 850
        self.assertAlmostEqual(result, 850.0, places=0)

    def test_sharp_feature_smoothed(self):
        """A sharp dip at center should be smoothed by band integration."""
        wls = np.arange(760, 781, 0.5)
        values = np.ones_like(wls)
        # Sharp O2-A absorption dip at 766nm only
        values[(wls >= 765) & (wls <= 767)] = 0.1

        # Center wavelength only would give ~0.1
        center_val = np.interp(766, wls, values)
        self.assertLess(center_val, 0.2)

        # Band-integrated with 6nm FWHM should average out the dip
        integrated = _band_integrate_result(766, 6.0, wls, values)
        self.assertGreater(integrated, center_val)


class TestFwhmFilePreference(unittest.TestCase):
    """Tests for find_coef_file() FWHM-specific file preference."""

    def setUp(self):
        """Create temporary COEFS directory with test files."""
        self.tmpdir = tempfile.mkdtemp(prefix='test_coefs_')
        self.atype_dir = os.path.join(self.tmpdir, 'CONTINENTAL')
        os.makedirs(self.atype_dir)

        # Create a generic file
        coef = SMACCoefficients()
        coef.a0T = 1.0  # Generic marker
        coef.to_file(os.path.join(self.atype_dir,
                                   'coef_766nm_CONTINENTAL.dat'))

        # Create an FWHM-specific file with different value
        coef.a0T = 0.98  # FWHM-specific marker
        coef.to_file(os.path.join(self.atype_dir,
                                   'coef_766nm_fwhm6nm_CONTINENTAL.dat'))

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_fwhm_file_preferred(self):
        """find_coef_file() should prefer FWHM-specific file when fwhm is set."""
        from radtran import find_coef_file
        path, wl = find_coef_file(766, 'continental', self.tmpdir, fwhm=6)
        self.assertIsNotNone(path)
        self.assertIn('fwhm6nm', path)

    def test_fallback_to_generic(self):
        """find_coef_file() should fall back to generic when no FWHM file exists."""
        from radtran import find_coef_file
        # Request fwhm=10 — no such file exists
        path, wl = find_coef_file(766, 'continental', self.tmpdir, fwhm=10)
        self.assertIsNotNone(path)
        self.assertNotIn('fwhm', path)

    def test_no_fwhm_uses_generic(self):
        """find_coef_file() should use generic when fwhm is None."""
        from radtran import find_coef_file
        path, wl = find_coef_file(766, 'continental', self.tmpdir, fwhm=None)
        self.assertIsNotNone(path)
        self.assertNotIn('fwhm', path)


class TestCommentLineSkipping(unittest.TestCase):
    """Tests for smac.coeff() skipping # comment lines."""

    def test_comment_lines_skipped(self):
        """Parser should handle # header lines in coefficient files."""
        from smac import coeff as SmacCoeff

        # Create a coefficient file with comment header
        coef = SMACCoefficients()
        coef.ah2o = -0.05
        coef.a0T = 0.97

        tmpdir = tempfile.mkdtemp(prefix='test_comments_')
        try:
            # Write file with comments
            plain_path = os.path.join(tmpdir, 'plain.dat')
            coef.to_file(plain_path)

            # Read plain file content and prepend comments
            with open(plain_path) as f:
                data_lines = f.read()

            comment_path = os.path.join(tmpdir, 'commented.dat')
            with open(comment_path, 'w') as f:
                f.write("# Band-integrated SMAC coefficients\n")
                f.write("# Center: 766nm, FWHM: 6nm, Type: CONTINENTAL\n")
                f.write(data_lines)

            # Both should parse identically
            c_plain = SmacCoeff(plain_path)
            c_comment = SmacCoeff(comment_path)

            self.assertAlmostEqual(c_plain.ah2o, c_comment.ah2o)
            self.assertAlmostEqual(c_plain.a0T, c_comment.a0T)
            self.assertAlmostEqual(c_plain.wo, c_comment.wo)
        finally:
            shutil.rmtree(tmpdir)


if __name__ == '__main__':
    unittest.main()
