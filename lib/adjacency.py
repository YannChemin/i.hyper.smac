#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Adjacency effect correction following Vermote et al. (1997).

In heterogeneous scenes, photons scattered from neighbouring pixels
contribute through the diffuse transmittance component.  Dark pixels
next to bright ones appear too bright; bright pixels next to dark ones
appear too dark.  This module implements a simple first-order correction
based on the environmental reflectance computed via a spatial average.
"""

import numpy as np
from scipy.ndimage import uniform_filter


def compute_environmental_reflectance(reflectance, pixel_size, psf_radius_km=1.0):
    """Compute spatially-averaged reflectance over a neighbourhood.

    Uses a box filter (uniform_filter) whose size matches the point
    spread function radius of atmospheric scattering.

    Args:
        reflectance: 2D array [rows, cols] of BOA reflectance.
        pixel_size: Pixel size in metres (assumes square pixels).
        psf_radius_km: Radius of the environmental PSF in km (default 1.0).

    Returns:
        2D array of environmental reflectance, same shape as input.
    """
    reflectance = np.asarray(reflectance, dtype=float)
    psf_radius_m = psf_radius_km * 1000.0
    # Filter size in pixels (diameter), minimum 3
    filter_size = max(3, int(round(2 * psf_radius_m / pixel_size)))
    # Make odd
    if filter_size % 2 == 0:
        filter_size += 1

    # Handle NaN by filling temporarily, then restoring
    nan_mask = np.isnan(reflectance)
    filled = reflectance.copy()
    if np.any(nan_mask):
        filled[nan_mask] = 0.0
        weight = (~nan_mask).astype(float)
        sum_vals = uniform_filter(filled, size=filter_size, mode='nearest')
        sum_wts = uniform_filter(weight, size=filter_size, mode='nearest')
        with np.errstate(divide='ignore', invalid='ignore'):
            r_env = sum_vals / np.maximum(sum_wts, 1e-10)
        r_env[nan_mask] = np.nan
    else:
        r_env = uniform_filter(filled, size=filter_size, mode='nearest')

    return r_env


def compute_direct_transmittance(wavelength, aod550, pressure, sza, vza):
    """Beer-Lambert direct transmittance (two-way).

    T_dir = exp(-tau / cos(sza)) * exp(-tau / cos(vza))

    where tau = tau_rayleigh + tau_aerosol (scattering only, no gas).

    Args:
        wavelength: Band centre wavelength in nm.
        aod550: AOD at 550 nm (scalar or 2D array).
        pressure: Atmospheric pressure in hPa (scalar or 2D array).
        sza: Solar zenith angle in degrees (scalar).
        vza: View zenith angle in degrees (scalar).

    Returns:
        Direct transmittance T_dir (same shape as aod550).
    """
    aod550 = np.asarray(aod550, dtype=float)
    pressure = np.asarray(pressure, dtype=float)

    wl_um = wavelength / 1000.0
    tau_r = 0.008569 * wl_um**(-4) * (1 + 0.0113 * wl_um**(-2))
    tau_r = tau_r * pressure / 1013.25
    tau_a = aod550 * (wavelength / 550.0)**(-1.3)
    tau = tau_r + tau_a

    us = np.cos(np.radians(sza))
    uv = np.cos(np.radians(vza))

    return np.exp(-tau / us) * np.exp(-tau / uv)


def adjacency_correction(r_boa, T_scat, s, wavelength, aod550, pressure,
                          sza, vza, pixel_size, psf_radius_km=1.0):
    """Vermote et al. (1997) adjacency effect correction.

    The correction replaces the assumption that all surface photons
    originate from the observed pixel.  The diffuse transmittance
    component picks up signal from the environmental reflectance
    (spatial average) rather than the pixel reflectance:

        T_diff = clip(T_scat - T_dir, 0, T_scat)
        r_env  = spatial_average(r_boa)
        correction = T_diff * s * (r_boa - r_env) / (1 - s*r_env + eps)
        r_adj  = r_boa + correction

    Args:
        r_boa: 2D array [rows, cols] of BOA reflectance.
        T_scat: Two-way total scattering transmittance (scalar or 2D).
        s: Spherical albedo (scalar or 2D).
        wavelength: Band centre wavelength in nm.
        aod550: AOD at 550 nm (scalar or 2D array).
        pressure: Atmospheric pressure in hPa (scalar or 2D array).
        sza: Solar zenith angle in degrees.
        vza: View zenith angle in degrees.
        pixel_size: Pixel size in metres.
        psf_radius_km: Environmental PSF radius in km.

    Returns:
        Corrected 2D reflectance array.
    """
    r_boa = np.asarray(r_boa, dtype=float)
    T_scat = np.asarray(T_scat, dtype=float)
    s = np.asarray(s, dtype=float)

    T_dir = compute_direct_transmittance(wavelength, aod550, pressure, sza, vza)
    T_diff = np.clip(T_scat - T_dir, 0.0, T_scat)

    r_env = compute_environmental_reflectance(r_boa, pixel_size, psf_radius_km)

    eps = 1e-10
    correction = T_diff * s * (r_boa - r_env) / (1.0 - s * r_env + eps)

    r_adj = r_boa + correction

    # Preserve NaN
    r_adj[np.isnan(r_boa)] = np.nan

    return r_adj
