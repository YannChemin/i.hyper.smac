#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Spectral polishing: moving-median outlier detection and quality weighting.

Detects and optionally replaces spectral outlier bands that arise from
residual errors at sharp gas absorption edges (O2-A, H2O 940nm, etc.).
"""

import numpy as np


def compute_quality_weights(wavelengths, gas_transmission, threshold=0.10):
    """Per-band quality weight in [0, 1] based on gas transmission.

    Low-transmission bands and bands near sharp absorption edges receive
    lower weights.  The formula is:

        quality = clip(tg / threshold, 0, 1) * exp(-10 * |d(tg)/d(wl)|)

    Args:
        wavelengths: 1D array of band centre wavelengths (nm), shape [n_bands].
        gas_transmission: 1D array of gas-only transmission per band.
        threshold: Transmission below which quality drops to zero.

    Returns:
        1D array of quality weights in [0, 1], shape [n_bands].
    """
    wavelengths = np.asarray(wavelengths, dtype=float)
    tg = np.asarray(gas_transmission, dtype=float)

    # Base quality from absolute transmission
    q_abs = np.clip(tg / threshold, 0.0, 1.0)

    # Edge penalty from spectral gradient of tg
    dtg = np.gradient(tg, wavelengths)
    q_edge = np.exp(-10.0 * np.abs(dtg))

    return q_abs * q_edge


def spectral_polish(reflectance, wavelengths, quality_weights=None,
                    window=7, mad_threshold=3.0, replace=True):
    """Moving-median outlier detection with MAD-based flagging.

    For each spectral position, computes the local median and MAD
    (median absolute deviation) over a sliding window.  Bands whose
    deviation exceeds ``mad_threshold * MAD`` are flagged.

    Args:
        reflectance: Array of shape [n_bands] (single pixel) or
            [n_pixels, n_bands] (batch).
        wavelengths: 1D array of band centre wavelengths, shape [n_bands].
        quality_weights: Optional 1D quality weights [n_bands].  When
            provided, bands with quality < 0.5 are automatically flagged.
        window: Half-window size for the moving median (in bands).
        mad_threshold: Number of MADs for outlier detection.
        replace: If True, replace flagged values with interpolated ones.

    Returns:
        Tuple (polished_reflectance, flags) where *flags* is a boolean
        array of the same shape (True = flagged as outlier).
    """
    reflectance = np.asarray(reflectance, dtype=float)
    wavelengths = np.asarray(wavelengths, dtype=float)
    single_pixel = reflectance.ndim == 1

    if single_pixel:
        reflectance = reflectance[np.newaxis, :]  # [1, n_bands]

    n_pixels, n_bands = reflectance.shape
    flags = np.zeros_like(reflectance, dtype=bool)

    # Mark NaN bands
    nan_mask = np.isnan(reflectance)

    for j in range(n_bands):
        lo = max(0, j - window)
        hi = min(n_bands, j + window + 1)
        local = reflectance[:, lo:hi].copy()

        # Replace NaN with inf so they don't affect median
        local[np.isnan(local)] = np.inf

        med = np.nanmedian(
            np.where(np.isinf(local), np.nan, local), axis=1
        )  # [n_pixels]
        mad = np.nanmedian(
            np.abs(np.where(np.isinf(local), np.nan, local)
                   - med[:, np.newaxis]),
            axis=1
        )
        mad = np.maximum(mad, 1e-8)  # avoid division by zero

        deviation = np.abs(reflectance[:, j] - med)
        flags[:, j] = deviation > mad_threshold * mad

    # Don't flag NaN bands (they are preserved separately)
    flags[nan_mask] = False

    # Flag low-quality bands
    if quality_weights is not None:
        qw = np.asarray(quality_weights, dtype=float)
        low_quality = qw < 0.5
        flags[:, low_quality] = True

    result = reflectance.copy()
    if replace:
        result = interpolate_flagged(result, flags, wavelengths)

    # Restore NaN
    result[nan_mask] = np.nan

    if single_pixel:
        return result[0], flags[0]
    return result, flags


def interpolate_flagged(reflectance, flags, wavelengths):
    """Replace flagged bands with linear interpolation from unflagged neighbours.

    Args:
        reflectance: Array [n_pixels, n_bands].
        flags: Boolean array [n_pixels, n_bands] (True = flagged).
        wavelengths: 1D array [n_bands].

    Returns:
        Array with flagged values replaced by interpolation.
    """
    reflectance = reflectance.copy()
    wavelengths = np.asarray(wavelengths, dtype=float)
    n_pixels, n_bands = reflectance.shape

    for i in range(n_pixels):
        flagged = flags[i]
        good = ~flagged & ~np.isnan(reflectance[i])
        if np.sum(good) < 2:
            continue  # not enough anchors for interpolation
        if not np.any(flagged):
            continue

        reflectance[i, flagged] = np.interp(
            wavelengths[flagged],
            wavelengths[good],
            reflectance[i, good],
        )

    return reflectance
