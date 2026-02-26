#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Superpixel atmospheric retrieval and spatial interpolation.

Retrieves atmospheric parameters (AOD, H2O) on superpixel means for higher
SNR, then spatially interpolates to every pixel.  Based on the ISOFIT
approach of using spatial structure to reduce noise in atmosphere estimates.
"""

import numpy as np


def segment_image(radiance_band, n_segments=200):
    """SLIC superpixel segmentation on a single reference band (e.g., 860nm NIR).

    Args:
        radiance_band: 2D array [rows, cols] of radiance or reflectance.
        n_segments: Approximate number of superpixels (default 200).

    Returns:
        Label array [rows, cols] with integer superpixel IDs.
    """
    from skimage.segmentation import slic

    # Normalize to [0, 1] for SLIC
    band = np.asarray(radiance_band, dtype=float)
    nan_mask = np.isnan(band)
    band_filled = band.copy()
    if np.any(nan_mask):
        band_filled[nan_mask] = np.nanmedian(band)

    vmin, vmax = np.nanmin(band_filled), np.nanmax(band_filled)
    if vmax > vmin:
        band_norm = (band_filled - vmin) / (vmax - vmin)
    else:
        band_norm = np.zeros_like(band_filled)

    # SLIC expects [rows, cols, channels] or [rows, cols] for single channel
    # Use channel_axis=None for grayscale
    labels = slic(band_norm, n_segments=n_segments, compactness=10.0,
                  start_label=0, channel_axis=None)

    return labels


def superpixel_means(data_2d, labels):
    """Compute mean of data_2d per superpixel label.

    Args:
        data_2d: 2D array [rows, cols] of values.
        labels: 2D array [rows, cols] of integer superpixel labels.

    Returns:
        Dict {label: mean_value}, excluding NaN pixels from averaging.
    """
    data_2d = np.asarray(data_2d, dtype=float)
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)

    result = {}
    for lbl in unique_labels:
        mask = labels == lbl
        values = data_2d[mask]
        finite = values[np.isfinite(values)]
        if len(finite) > 0:
            result[lbl] = np.mean(finite)
        else:
            result[lbl] = np.nan

    return result


def interpolate_superpixel_field(labels, superpixel_values, smoothing_sigma=2.0):
    """Convert per-superpixel values to a smooth per-pixel field.

    1. Map superpixel_values to pixel grid via labels
    2. Apply Gaussian smoothing for spatial continuity

    Args:
        labels: 2D array [rows, cols] of superpixel labels.
        superpixel_values: Dict {label: value} from superpixel_means().
        smoothing_sigma: Gaussian smoothing sigma in pixels (default 2.0).

    Returns:
        2D array [rows, cols] of spatially smooth values.
    """
    from scipy.ndimage import gaussian_filter

    labels = np.asarray(labels)
    field = np.full(labels.shape, np.nan, dtype=float)

    for lbl, val in superpixel_values.items():
        field[labels == lbl] = val

    # Fill NaN with nearest valid before smoothing
    nan_mask = np.isnan(field)
    if np.any(nan_mask) and not np.all(nan_mask):
        from scipy.ndimage import distance_transform_edt
        _, nearest_idx = distance_transform_edt(nan_mask, return_distances=True,
                                                 return_indices=True)
        field[nan_mask] = field[tuple(nearest_idx[:, nan_mask])]

    if smoothing_sigma > 0:
        field = gaussian_filter(field, sigma=smoothing_sigma)

    return field
