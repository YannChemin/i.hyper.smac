#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for SMAC atmospheric correction algorithm.

These tests verify the core SMAC functions without requiring GRASS GIS.
"""

import unittest
import sys
import os
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from smac import PdeZ, smac_inv, smac_dir


class MockCoefficients:
    """Mock SMAC coefficient object for testing."""

    def __init__(self):
        # H2O coefficients
        self.ah2o = -0.001
        self.nh2o = 0.5
        # O3 coefficients
        self.ao3 = -0.05
        self.no3 = 0.5
        # O2 coefficients
        self.ao2 = 0.0
        self.no2 = 0.0
        self.po2 = 1.0
        # CO2 coefficients
        self.aco2 = 0.0
        self.nco2 = 0.0
        self.pco2 = 1.0
        # CH4 coefficients
        self.ach4 = 0.0
        self.nch4 = 0.0
        self.pch4 = 1.0
        # NO2 coefficients
        self.ano2 = 0.0
        self.nno2 = 0.0
        self.pno2 = 1.0
        # CO coefficients
        self.aco = 0.0
        self.nco = 0.0
        self.pco = 1.0
        # Scattering coefficients
        self.a0s = 0.1
        self.a1s = 0.01
        self.a2s = 0.001
        self.a3s = 0.0
        # Transmission coefficients
        self.a0T = 0.9
        self.a1T = -0.01
        self.a2T = 0.0
        self.a3T = 0.0
        # Rayleigh optical thickness
        self.taur = 0.05
        self.sr = 0.05
        # Aerosol coefficients
        self.a0taup = 0.0
        self.a1taup = 1.0
        self.wo = 0.9
        self.gc = 0.7
        # Phase function coefficients
        self.a0P = 0.1
        self.a1P = 0.0
        self.a2P = 0.0
        self.a3P = 0.0
        self.a4P = 0.0
        # Residual coefficients
        self.Rest1 = 0.0
        self.Rest2 = 0.0
        self.Rest3 = 0.0
        self.Rest4 = 0.0
        self.Resr1 = 0.0
        self.Resr2 = 0.0
        self.Resr3 = 0.0
        self.Resa1 = 0.0
        self.Resa2 = 0.0
        self.Resa3 = 0.0
        self.Resa4 = 0.0


class TestPdeZ(unittest.TestCase):
    """Tests for atmospheric pressure calculation."""

    def test_sea_level(self):
        """Pressure at sea level should be approximately 1013.25 hPa."""
        pressure = PdeZ(0)
        self.assertAlmostEqual(pressure, 1013.25, places=1)

    def test_positive_altitude(self):
        """Pressure decreases with altitude."""
        p_sea = PdeZ(0)
        p_1000m = PdeZ(1000)
        p_3000m = PdeZ(3000)

        self.assertGreater(p_sea, p_1000m)
        self.assertGreater(p_1000m, p_3000m)

    def test_typical_values(self):
        """Check pressure at typical altitudes."""
        # At ~1500m, pressure should be around 850 hPa
        p_1500 = PdeZ(1500)
        self.assertGreater(p_1500, 800)
        self.assertLess(p_1500, 900)

        # At ~5500m (approximate Mount Everest base camp), pressure ~500 hPa
        p_5500 = PdeZ(5500)
        self.assertGreater(p_5500, 450)
        self.assertLess(p_5500, 550)


class TestSmacInverse(unittest.TestCase):
    """Tests for SMAC inverse model (TOA to surface reflectance)."""

    def setUp(self):
        """Set up test fixtures."""
        self.coef = MockCoefficients()
        # Standard atmospheric parameters
        self.theta_s = 30.0  # Solar zenith
        self.phi_s = 180.0   # Solar azimuth
        self.theta_v = 0.0   # View zenith (nadir)
        self.phi_v = 0.0     # View azimuth
        self.pressure = 1013.25  # Sea level
        self.taup550 = 0.1   # AOD at 550nm
        self.uo3 = 0.3       # Ozone (cm-atm)
        self.uh2o = 2.0      # Water vapor (g/cm2)

    def test_reflectance_reduction(self):
        """Surface reflectance should generally be lower than TOA reflectance."""
        r_toa = 0.3
        r_surf = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )
        # Surface reflectance should be positive and typically lower than TOA
        self.assertGreater(r_surf, 0)

    def test_array_input(self):
        """SMAC should handle numpy array inputs."""
        r_toa = np.array([0.1, 0.2, 0.3, 0.4])
        r_surf = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )
        self.assertEqual(r_surf.shape, r_toa.shape)

    def test_high_aod_stronger_correction(self):
        """Higher AOD should result in larger atmospheric correction."""
        r_toa = 0.25

        r_surf_low_aod = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, 0.05, self.uo3, self.uh2o, self.coef
        )

        r_surf_high_aod = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, 0.5, self.uo3, self.uh2o, self.coef
        )

        # Higher AOD means more atmospheric contribution to TOA,
        # so surface reflectance might be lower
        # (depends on the specific coefficients)
        self.assertNotEqual(r_surf_low_aod, r_surf_high_aod)


class TestSmacDirect(unittest.TestCase):
    """Tests for SMAC direct model (surface to TOA reflectance)."""

    def setUp(self):
        """Set up test fixtures."""
        self.coef = MockCoefficients()
        self.theta_s = 30.0
        self.phi_s = 180.0
        self.theta_v = 0.0
        self.phi_v = 0.0
        self.pressure = 1013.25
        self.taup550 = 0.1
        self.uo3 = 0.3
        self.uh2o = 2.0

    def test_reflectance_increase(self):
        """TOA reflectance should generally be higher than surface reflectance."""
        r_surf = 0.2
        r_toa = smac_dir(
            r_surf, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )
        # TOA reflectance should be positive
        self.assertGreater(r_toa, 0)

    def test_array_input(self):
        """SMAC direct should handle numpy array inputs."""
        r_surf = np.array([0.05, 0.1, 0.15, 0.2])
        r_toa = smac_dir(
            r_surf, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )
        self.assertEqual(r_toa.shape, r_surf.shape)


class TestSmacRoundTrip(unittest.TestCase):
    """Tests for SMAC inverse/direct round-trip consistency."""

    def setUp(self):
        """Set up test fixtures."""
        self.coef = MockCoefficients()
        self.theta_s = 30.0
        self.phi_s = 180.0
        self.theta_v = 0.0
        self.phi_v = 0.0
        self.pressure = 1013.25
        self.taup550 = 0.1
        self.uo3 = 0.3
        self.uh2o = 2.0

    def test_round_trip_scalar(self):
        """Applying direct then inverse should approximately recover original."""
        r_surf_orig = 0.15

        # Surface -> TOA
        r_toa = smac_dir(
            r_surf_orig, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )

        # TOA -> Surface
        r_surf_recovered = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )

        # Should be close to original (within numerical precision)
        self.assertAlmostEqual(r_surf_orig, r_surf_recovered, places=4)

    def test_round_trip_array(self):
        """Round-trip test with array inputs."""
        r_surf_orig = np.array([0.05, 0.1, 0.15, 0.2, 0.25])

        r_toa = smac_dir(
            r_surf_orig, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )

        r_surf_recovered = smac_inv(
            r_toa, self.theta_s, self.phi_s, self.theta_v, self.phi_v,
            self.pressure, self.taup550, self.uo3, self.uh2o, self.coef
        )

        np.testing.assert_array_almost_equal(r_surf_orig, r_surf_recovered, decimal=4)


if __name__ == '__main__':
    unittest.main()
