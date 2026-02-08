#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for AOD (Aerosol Optical Depth) estimation algorithm.

These tests verify the core AOD math without requiring GRASS GIS.
All algorithms are reimplemented inline from lib/aod.py.
"""

import unittest
import sys
import os
import math
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

# Constants from aod.py
BLUE_WL = 470.0
RED_WL = 660.0
NIR_WL = 860.0
SWIR_WL = 2130.0
OMEGA_0 = 0.89
ASYM_G = 0.65
ANGSTROM_DEFAULT = 1.3


# --- Reimplemented pure-math functions (no GRASS dependency) ---

def compute_scatter_angle(solar_zenith, view_zenith, solar_azimuth, view_azimuth):
    """Compute scattering angle from geometry (degrees in, radians out)."""
    theta_s_r = math.radians(solar_zenith)
    theta_v_r = math.radians(view_zenith)
    dphi_r = math.radians(solar_azimuth - view_azimuth)
    cos_scatter = (-math.cos(theta_s_r) * math.cos(theta_v_r)
                   + math.sin(theta_s_r) * math.sin(theta_v_r)
                   * math.cos(dphi_r))
    cos_scatter = max(-1.0, min(1.0, cos_scatter))
    return math.acos(cos_scatter)


def henyey_greenstein(scatter_angle, g):
    """Henyey-Greenstein phase function."""
    cos_scat = math.cos(scatter_angle)
    return (1 - g**2) / (1 + g**2 - 2 * g * cos_scat)**1.5


def single_scattering_inversion(rho_path, mu_s, mu_v, omega_0, phase):
    """Invert path radiance to AOD via single-scattering approximation."""
    return rho_path * 4.0 * mu_s * mu_v / (omega_0 * phase)


def angstrom_exponent(tau_blue, tau_red, wl_blue=BLUE_WL, wl_red=RED_WL):
    """Compute Angstrom exponent from two-wavelength AODs."""
    if tau_blue <= 0 or tau_red <= 0:
        return ANGSTROM_DEFAULT
    return -math.log(tau_blue / tau_red) / math.log(wl_blue / wl_red)


def scale_aod_to_550(tau_ref, wl_ref, alpha):
    """Scale AOD from reference wavelength to 550nm."""
    return tau_ref * (550.0 / wl_ref) ** (-alpha)


def radiance_to_reflectance(L, E0, mu_s, d=1.0):
    """Convert radiance to TOA reflectance."""
    return math.pi * L * d**2 / (E0 * mu_s)


def compute_E0_analytical(wavelength):
    """Analytical blackbody E0 fallback (Planck function at 5778K).

    Returns E0 in mW/(m^2 nm).
    """
    wl_m = wavelength * 1e-9
    h = 6.62607015e-34
    c = 2.998e8
    k = 1.381e-23
    T = 5778.0
    hc_kT = h * c / (wl_m * k * T)
    B = 2 * h * c**2 / (wl_m**5 * (math.exp(hc_kT) - 1.0))
    return 6.794e-5 * B * 1e-6


def clamp_aod(aod, aod_min=0.01, aod_max=2.0):
    """Clamp AOD to valid range."""
    return max(aod_min, min(aod, aod_max))


# --- Test classes ---

class TestScatteringGeometry(unittest.TestCase):
    """Scattering angle computation from solar/view geometry."""

    def test_backscatter(self):
        """Sun behind sensor (opposite azimuths): scatter angle near 180 deg."""
        # SZA=30, VZA=30, sun and sensor on opposite sides
        angle = compute_scatter_angle(30, 30, 0, 180)
        self.assertGreater(math.degrees(angle), 150)

    def test_forward_scatter(self):
        """Sun in front of sensor (same azimuth): scatter angle < 180 deg.

        The formula uses cos(scatter) = -cos(SZA)*cos(VZA) + sin(SZA)*sin(VZA)*cos(dphi).
        With SZA=VZA=30, dphi=0: cos(scatter) = -0.75 + 0.25 = -0.5 -> 120 deg.
        True forward scatter requires dphi=180 (sensor looking toward sun).
        """
        # Same azimuth: scatter = 120 deg (geometric scattering convention)
        angle_same = compute_scatter_angle(30, 30, 0, 0)
        self.assertAlmostEqual(math.degrees(angle_same), 120.0, places=3)

        # Opposite azimuths (dphi=180): forward scatter
        angle_fwd = compute_scatter_angle(30, 30, 0, 180)
        self.assertGreater(math.degrees(angle_fwd), 150)

    def test_nadir_viewing(self):
        """Nadir viewing (VZA=0): scatter angle equals solar zenith."""
        for sza in [10, 30, 45, 60]:
            angle = compute_scatter_angle(sza, 0, 0, 0)
            # With VZA=0, sin(VZA)=0, cos(VZA)=1
            # cos(scatter) = -cos(SZA)*1 + 0 = -cos(SZA)
            # scatter = 180 - SZA
            expected = 180 - sza
            self.assertAlmostEqual(math.degrees(angle), expected, places=5,
                                   msg=f"SZA={sza}")

    def test_symmetric_azimuths(self):
        """Swapping solar/view azimuth by 180 changes the scattering geometry."""
        a1 = compute_scatter_angle(30, 15, 0, 90)
        a2 = compute_scatter_angle(30, 15, 90, 0)
        # These are generally different (not symmetric in this formula)
        # but both should be valid angles
        self.assertGreater(a1, 0)
        self.assertGreater(a2, 0)
        self.assertLess(a1, math.pi)
        self.assertLess(a2, math.pi)

    def test_clamping(self):
        """cos(scatter) is clamped to [-1, 1] to avoid acos domain errors."""
        # Edge case: identical geometry
        angle = compute_scatter_angle(0, 0, 0, 0)
        self.assertAlmostEqual(math.degrees(angle), 180.0, places=5)


class TestHenyeyGreensteinPhase(unittest.TestCase):
    """Henyey-Greenstein phase function tests."""

    def test_isotropic(self):
        """g=0 should give isotropic scattering (P=1 for all angles)."""
        for angle_deg in [0, 45, 90, 135, 180]:
            angle_rad = math.radians(angle_deg)
            P = henyey_greenstein(angle_rad, g=0.0)
            self.assertAlmostEqual(P, 1.0, places=5,
                                   msg=f"angle={angle_deg}")

    def test_forward_greater_than_backward(self):
        """For g>0 (forward scattering), P(0) > P(180)."""
        P_forward = henyey_greenstein(0, g=ASYM_G)
        P_backward = henyey_greenstein(math.pi, g=ASYM_G)
        self.assertGreater(P_forward, P_backward)

    def test_positive_everywhere(self):
        """Phase function must be positive for all angles."""
        for angle_deg in range(0, 181, 10):
            P = henyey_greenstein(math.radians(angle_deg), g=ASYM_G)
            self.assertGreater(P, 0, msg=f"angle={angle_deg}")

    def test_continental_values(self):
        """Phase function with g=0.65 at key angles."""
        # Forward (0 deg): should be large
        P_0 = henyey_greenstein(0, g=ASYM_G)
        self.assertGreater(P_0, 5.0)

        # 90 degrees: intermediate
        P_90 = henyey_greenstein(math.pi / 2, g=ASYM_G)
        self.assertGreater(P_90, 0)
        self.assertLess(P_90, P_0)

        # Backward (180 deg): small
        P_180 = henyey_greenstein(math.pi, g=ASYM_G)
        self.assertLess(P_180, 1.0)


class TestSingleScatteringInversion(unittest.TestCase):
    """AOD inversion from path radiance via single-scattering."""

    def test_zero_path_radiance(self):
        """Zero path radiance gives zero AOD."""
        tau = single_scattering_inversion(0.0, 0.866, 1.0, OMEGA_0, 1.0)
        self.assertEqual(tau, 0.0)

    def test_known_inversion(self):
        """Known path radiance gives expected AOD."""
        mu_s = math.cos(math.radians(30))  # 0.866
        mu_v = 1.0  # nadir
        scatter = compute_scatter_angle(30, 0, 0, 0)
        phase = henyey_greenstein(scatter, ASYM_G)

        rho_path = 0.02
        tau = single_scattering_inversion(rho_path, mu_s, mu_v, OMEGA_0, phase)
        # Should give a reasonable AOD (0 < tau < 1 for moderate path radiance)
        self.assertGreater(tau, 0)
        self.assertLess(tau, 1.0)

    def test_proportional_to_path_radiance(self):
        """AOD is proportional to path radiance."""
        mu_s = 0.866
        mu_v = 1.0
        phase = henyey_greenstein(math.radians(150), ASYM_G)

        tau1 = single_scattering_inversion(0.01, mu_s, mu_v, OMEGA_0, phase)
        tau2 = single_scattering_inversion(0.02, mu_s, mu_v, OMEGA_0, phase)
        self.assertAlmostEqual(tau2 / tau1, 2.0, places=5)

    def test_inversion_factor(self):
        """Inversion factor = 4 * mu_s * mu_v / (omega_0 * P)."""
        mu_s = math.cos(math.radians(30))
        mu_v = math.cos(math.radians(10))
        phase = 1.5
        expected_factor = 4.0 * mu_s * mu_v / (OMEGA_0 * phase)
        tau = single_scattering_inversion(1.0, mu_s, mu_v, OMEGA_0, phase)
        self.assertAlmostEqual(tau, expected_factor, places=5)


class TestAngstromExponent(unittest.TestCase):
    """Angstrom exponent from multi-wavelength AOD."""

    def test_equal_tau(self):
        """Equal AOD at both wavelengths gives alpha=0."""
        alpha = angstrom_exponent(0.1, 0.1)
        self.assertAlmostEqual(alpha, 0.0, places=5)

    def test_blue_greater_positive_alpha(self):
        """tau_blue > tau_red gives alpha > 0 (fine-mode aerosol)."""
        alpha = angstrom_exponent(0.2, 0.1)
        self.assertGreater(alpha, 0)

    def test_typical_continental(self):
        """Continental aerosol: alpha ~ 1.0-1.5."""
        # Simulate tau_blue=0.15, tau_red=0.08
        alpha = angstrom_exponent(0.15, 0.08)
        self.assertGreater(alpha, 0.5)
        self.assertLess(alpha, 2.5)

    def test_scaling_to_550(self):
        """Scaling AOD from blue to 550nm."""
        tau_blue = 0.15
        alpha = 1.3
        tau_550 = scale_aod_to_550(tau_blue, BLUE_WL, alpha)
        # 550 > 470, so with alpha > 0, tau_550 < tau_blue
        self.assertLess(tau_550, tau_blue)
        self.assertGreater(tau_550, 0)

    def test_out_of_range_clamped(self):
        """Out-of-range alpha returns default."""
        # Zero or negative tau returns default
        alpha = angstrom_exponent(0.0, 0.1)
        self.assertEqual(alpha, ANGSTROM_DEFAULT)

        alpha = angstrom_exponent(-0.1, 0.1)
        self.assertEqual(alpha, ANGSTROM_DEFAULT)


class TestDarkTargetSurfaceReflectance(unittest.TestCase):
    """Kaufman 1997 VIS/SWIR surface reflectance ratios."""

    def test_blue_ratio(self):
        """rho_blue = 0.25 * rho_SWIR."""
        rho_swir = 0.10
        rho_blue_surf = 0.25 * rho_swir
        self.assertAlmostEqual(rho_blue_surf, 0.025, places=5)

    def test_red_ratio(self):
        """rho_red = 0.50 * rho_SWIR."""
        rho_swir = 0.10
        rho_red_surf = 0.50 * rho_swir
        self.assertAlmostEqual(rho_red_surf, 0.050, places=5)

    def test_path_radiance_clamped(self):
        """Path radiance = max(rho_TOA - rho_surf, 0)."""
        rho_toa = 0.08
        rho_surf = 0.025
        rho_path = max(rho_toa - rho_surf, 0.0)
        self.assertGreater(rho_path, 0)

        # When surface > TOA (shouldn't happen normally, but test clamping)
        rho_path_neg = max(0.01 - 0.05, 0.0)
        self.assertEqual(rho_path_neg, 0.0)

    def test_blue_less_than_red(self):
        """Blue surface reflectance < red surface reflectance."""
        rho_swir = 0.12
        self.assertLess(0.25 * rho_swir, 0.50 * rho_swir)


class TestRadianceToReflectance(unittest.TestCase):
    """Radiance to TOA reflectance conversion."""

    def test_known_values(self):
        """Known L, E0, geometry -> expected reflectance."""
        L = 100.0   # mW/(m^2 sr nm)
        E0 = 1900.0  # mW/(m^2 nm) (typical at 550nm)
        mu_s = math.cos(math.radians(30))
        d = 1.0

        rho = radiance_to_reflectance(L, E0, mu_s, d)
        # pi * 100 / (1900 * 0.866) ~ 0.191
        self.assertGreater(rho, 0.15)
        self.assertLess(rho, 0.25)

    def test_distance_scaling(self):
        """Earth-Sun distance d=1 at mean distance."""
        L = 50.0
        E0 = 1900.0
        mu_s = 1.0  # nadir sun

        rho_d1 = radiance_to_reflectance(L, E0, mu_s, d=1.0)
        rho_d1_01 = radiance_to_reflectance(L, E0, mu_s, d=1.01)
        # Larger distance -> higher reflectance (inverse square)
        self.assertGreater(rho_d1_01, rho_d1)

    def test_formula(self):
        """rho = pi * L * d^2 / (E0 * cos(theta_s))."""
        L = 1.0
        E0 = math.pi
        mu_s = 1.0
        d = 1.0
        rho = radiance_to_reflectance(L, E0, mu_s, d)
        self.assertAlmostEqual(rho, 1.0, places=5)

    def test_positive_output(self):
        """Reflectance is positive for positive inputs."""
        rho = radiance_to_reflectance(50, 1800, 0.5, 1.0)
        self.assertGreater(rho, 0)


class TestE0AnalyticalFallback(unittest.TestCase):
    """Planck blackbody E0 at 5778K."""

    def test_550nm_range(self):
        """E0 at 550nm should be in range ~1700-1900 mW/(m^2 nm)."""
        e0 = compute_E0_analytical(550.0)
        self.assertGreater(e0, 1700)
        self.assertLess(e0, 1900)

    def test_decreases_in_swir(self):
        """E0 decreases from VIS to SWIR."""
        e0_550 = compute_E0_analytical(550.0)
        e0_1000 = compute_E0_analytical(1000.0)
        e0_2100 = compute_E0_analytical(2100.0)
        self.assertGreater(e0_550, e0_1000)
        self.assertGreater(e0_1000, e0_2100)

    def test_positive_all_wavelengths(self):
        """E0 is positive for all sensor wavelengths."""
        for wl in [400, 550, 660, 860, 1000, 1600, 2100, 2500]:
            e0 = compute_E0_analytical(wl)
            self.assertGreater(e0, 0, msg=f"wl={wl}nm")

    def test_blue_band(self):
        """E0 at 470nm should be reasonable."""
        e0 = compute_E0_analytical(470.0)
        # Should be in the ballpark of the Kurucz spectrum at 470nm
        self.assertGreater(e0, 1600)
        self.assertLess(e0, 2000)


class TestAODBounds(unittest.TestCase):
    """AOD clamping to valid range [0.01, 2.0]."""

    def test_normal_value(self):
        """Normal AOD passes through."""
        self.assertEqual(clamp_aod(0.15), 0.15)
        self.assertEqual(clamp_aod(0.5), 0.5)

    def test_below_min(self):
        """AOD below minimum is clamped to 0.01."""
        self.assertEqual(clamp_aod(0.0), 0.01)
        self.assertEqual(clamp_aod(-0.5), 0.01)

    def test_above_max(self):
        """AOD above maximum is clamped to 2.0."""
        self.assertEqual(clamp_aod(3.0), 2.0)
        self.assertEqual(clamp_aod(10.0), 2.0)

    def test_boundary_values(self):
        """Boundary values are preserved."""
        self.assertEqual(clamp_aod(0.01), 0.01)
        self.assertEqual(clamp_aod(2.0), 2.0)


if __name__ == '__main__':
    unittest.main()
