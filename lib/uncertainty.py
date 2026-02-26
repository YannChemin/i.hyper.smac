#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reflectance uncertainty propagation and model discrepancy.

Computes per-band, per-pixel reflectance uncertainty from three sources:
1. Instrument noise (NEDL)
2. AOD retrieval uncertainty
3. Systematic RT model errors (model discrepancy)

Based on the ISOFIT posterior uncertainty framework.
"""

import numpy as np


def estimate_nedl(radiance_band, percentile=5):
    """Estimate noise-equivalent delta radiance from dark pixels.

    Uses the standard deviation of radiance in the darkest pixels
    (below the given percentile) as a proxy for instrument noise.

    Args:
        radiance_band: 2D array of TOA radiance.
        percentile: Percentile threshold for "dark" pixels (default 5).

    Returns:
        Scalar NEDL estimate.
    """
    radiance = np.asarray(radiance_band, dtype=float)
    finite = radiance[np.isfinite(radiance) & (radiance > 0)]
    if len(finite) < 10:
        return 0.01  # fallback

    threshold = np.percentile(finite, percentile)
    dark = finite[finite <= threshold]
    if len(dark) < 5:
        return 0.01

    return float(np.std(dark))


def compute_reflectance_uncertainty(rad_toa, refl_boa, E0, d2, cos_sza,
                                     T_down, T_up, s, R_atm,
                                     nedl=None, aod_sigma=0.04,
                                     atm_lut=None, wavelength=None,
                                     aod=None, h2o=None):
    """Compute per-pixel reflectance uncertainty.

    Three sources propagated through the algebraic inversion:
    1. Instrument noise (NEDL -> sigma_rho_toa -> sigma_rfl)
    2. AOD uncertainty (delta_aod -> delta_R_atm, delta_T -> delta_rfl)
    3. H2O uncertainty is captured implicitly via superpixel smoothing

    Args:
        rad_toa: 2D array of TOA radiance.
        refl_boa: 2D array of retrieved BOA reflectance.
        E0: Exo-atmospheric irradiance (scalar).
        d2: Earth-Sun distance squared.
        cos_sza: Cosine of solar zenith angle.
        T_down: Downward transmittance (scalar or 2D).
        T_up: Upward transmittance (scalar or 2D).
        s: Spherical albedo (scalar or 2D).
        R_atm: Path reflectance (scalar or 2D).
        nedl: Noise-equivalent delta radiance. If None, estimated from dark pixels.
        aod_sigma: AOD uncertainty (default 0.04, ~20% of typical AOD).
        atm_lut: AtmosphericLUT instance for AOD perturbation.
        wavelength: Band wavelength (nm), needed for AOD perturbation.
        aod: Per-pixel or scalar AOD, needed for AOD perturbation.
        h2o: Per-pixel or scalar H2O, needed for LUT interpolation.

    Returns:
        2D array of sigma_rfl (reflectance uncertainty), same shape as refl_boa.
    """
    rad_toa = np.asarray(rad_toa, dtype=float)
    refl_boa = np.asarray(refl_boa, dtype=float)
    T_total = np.asarray(T_down, dtype=float) * np.asarray(T_up, dtype=float)

    # 1. Instrument noise contribution
    if nedl is None:
        nedl = estimate_nedl(rad_toa)

    # Convert NEDL to TOA reflectance noise
    sigma_rho_toa = np.pi * nedl * d2 / (E0 * cos_sza)

    # Propagate through inversion: first order
    sigma_rfl_noise = sigma_rho_toa / np.maximum(T_total, 1e-6)

    # 2. AOD uncertainty contribution
    sigma_rfl_aod = np.zeros_like(refl_boa)
    if atm_lut is not None and wavelength is not None and aod is not None:
        aod_arr = np.asarray(aod, dtype=float)
        aod_plus = aod_arr + aod_sigma
        aod_minus = np.maximum(aod_arr - aod_sigma, 0.001)

        R1, Td1, Tu1, s1 = atm_lut.interpolate(wavelength, aod_plus, h2o=h2o)
        R0, Td0, Tu0, s0 = atm_lut.interpolate(wavelength, aod_minus, h2o=h2o)

        # Invert at perturbed AOD
        refl_toa = np.asarray(R_atm, dtype=float) + T_total * refl_boa / (
            1.0 - np.asarray(s, dtype=float) * refl_boa + 1e-10)

        T1 = Td1 * Tu1
        y1 = (refl_toa - R1) / np.maximum(T1, 1e-6)
        rfl1 = y1 / (1.0 + s1 * y1 + 1e-10)

        T0 = Td0 * Tu0
        y0 = (refl_toa - R0) / np.maximum(T0, 1e-6)
        rfl0 = y0 / (1.0 + s0 * y0 + 1e-10)

        sigma_rfl_aod = np.abs(rfl1 - rfl0) / 2.0

    # Total uncertainty (RSS)
    sigma_total = np.sqrt(sigma_rfl_noise ** 2 + sigma_rfl_aod ** 2)

    # Preserve NaN
    sigma_total[np.isnan(refl_boa)] = np.nan

    return sigma_total


def compute_model_discrepancy(wavelengths):
    """Per-band model discrepancy (systematic RT error floor).

    Empirical characterization of where the libRadtran LUT is least accurate:
    - Gas absorption band edges (720, 760, 940, 1135nm): high discrepancy
    - SWIR (>1500nm): moderate discrepancy (aerosol model uncertainty)
    - VIS clear windows: low discrepancy

    Args:
        wavelengths: 1D array of band center wavelengths (nm).

    Returns:
        1D array of sigma_model per band (in reflectance units).
    """
    wavelengths = np.asarray(wavelengths, dtype=float)
    sigma = np.full_like(wavelengths, 0.005)  # baseline 0.5%

    # Higher discrepancy at gas absorption band edges
    for center in [720, 760, 940, 1135, 1380, 1850]:
        sigma += 0.02 * np.exp(-((wavelengths - center) / 30) ** 2)

    # Moderate increase in SWIR
    swir_mask = wavelengths > 1500
    sigma[swir_mask] += 0.01

    return sigma
