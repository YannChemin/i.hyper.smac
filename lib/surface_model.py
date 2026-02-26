#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Surface reflectance prior model (Gaussian mixture).

Provides a simple 3-component surface prior (vegetation, soil, water) for
MAP regularization of retrieved reflectance.  Based on the ISOFIT approach
of using surface knowledge to constrain the inversion in noisy bands.
"""

import numpy as np

# Built-in reference spectra at key wavelengths (nm -> reflectance)
# Derived from typical field/lab measurements.
VEGETATION_SPECTRUM = {
    400: 0.03, 450: 0.04, 500: 0.05, 550: 0.08, 600: 0.05,
    650: 0.04, 680: 0.03, 700: 0.08, 720: 0.20, 750: 0.40,
    800: 0.43, 850: 0.45, 900: 0.44, 1000: 0.42, 1100: 0.40,
    1200: 0.38, 1300: 0.35, 1400: 0.10, 1500: 0.30, 1600: 0.28,
    1700: 0.25, 1800: 0.10, 1900: 0.05, 2000: 0.20, 2100: 0.22,
    2200: 0.18, 2300: 0.12, 2400: 0.08, 2500: 0.05,
}

SOIL_SPECTRUM = {
    400: 0.04, 450: 0.05, 500: 0.07, 550: 0.10, 600: 0.13,
    650: 0.16, 680: 0.18, 700: 0.19, 720: 0.20, 750: 0.21,
    800: 0.22, 850: 0.24, 900: 0.25, 1000: 0.27, 1100: 0.28,
    1200: 0.30, 1300: 0.31, 1400: 0.15, 1500: 0.30, 1600: 0.32,
    1700: 0.33, 1800: 0.15, 1900: 0.10, 2000: 0.30, 2100: 0.31,
    2200: 0.28, 2300: 0.22, 2400: 0.18, 2500: 0.15,
}

WATER_SPECTRUM = {
    400: 0.05, 450: 0.06, 500: 0.06, 550: 0.05, 600: 0.04,
    650: 0.03, 680: 0.03, 700: 0.02, 720: 0.02, 750: 0.02,
    800: 0.01, 850: 0.01, 900: 0.01, 1000: 0.01, 1100: 0.005,
    1200: 0.005, 1300: 0.003, 1400: 0.002, 1500: 0.002, 1600: 0.001,
    1700: 0.001, 1800: 0.001, 1900: 0.001, 2000: 0.001, 2100: 0.001,
    2200: 0.001, 2300: 0.001, 2400: 0.001, 2500: 0.001,
}

# Per-component variance (squared standard deviation) around the mean.
# Larger variance = weaker prior constraint.  Units: reflectance^2.
VEGETATION_VARIANCE_SCALE = 0.15  # moderate variability
SOIL_VARIANCE_SCALE = 0.20        # high variability (many soil types)
WATER_VARIANCE_SCALE = 0.03       # low variability (water is dark)

REFERENCE_SPECTRA = [VEGETATION_SPECTRUM, SOIL_SPECTRUM, WATER_SPECTRUM]
VARIANCE_SCALES = [VEGETATION_VARIANCE_SCALE, SOIL_VARIANCE_SCALE,
                   WATER_VARIANCE_SCALE]


class SurfaceModel:
    """Simple 3-component surface prior for MAP regularization."""

    def __init__(self, wavelengths, n_components=3):
        """Initialize with sensor band wavelengths.

        Interpolates reference spectra to sensor wavelengths and computes
        per-component mean and diagonal covariance (variance per band).

        Args:
            wavelengths: 1D array of band center wavelengths (nm).
            n_components: Number of surface components (default 3).
        """
        self.wavelengths = np.asarray(wavelengths, dtype=float)
        self.n_bands = len(self.wavelengths)
        self.n_components = min(n_components, len(REFERENCE_SPECTRA))

        # Interpolate reference spectra to sensor wavelengths
        self.means = np.zeros((self.n_components, self.n_bands))
        self.variances = np.zeros((self.n_components, self.n_bands))

        for c in range(self.n_components):
            ref = REFERENCE_SPECTRA[c]
            ref_wl = np.array(sorted(ref.keys()), dtype=float)
            ref_val = np.array([ref[w] for w in sorted(ref.keys())], dtype=float)
            self.means[c] = np.interp(self.wavelengths, ref_wl, ref_val)
            # Variance = (scale * mean)^2, with a floor
            scale = VARIANCE_SCALES[c]
            self.variances[c] = np.maximum(
                (scale * self.means[c]) ** 2,
                0.001 ** 2  # minimum variance floor
            )

        # VNIR mask for classification (400-1000nm)
        self.vnir_mask = (self.wavelengths >= 400) & (self.wavelengths <= 1000)

    def classify_pixel(self, reflectance):
        """Return component index with minimum Euclidean distance.

        Uses VNIR bands only (400-1000nm) for robust classification,
        since these are least affected by atmospheric residuals.

        Args:
            reflectance: 1D array [n_bands] of reflectance values.

        Returns:
            Integer component index (0=vegetation, 1=soil, 2=water).
        """
        reflectance = np.asarray(reflectance, dtype=float)
        vnir = self.vnir_mask

        if not np.any(vnir) or np.all(np.isnan(reflectance[vnir])):
            return 1  # default to soil

        # Normalize by mean to reduce scale effects
        obs = reflectance[vnir]
        finite = np.isfinite(obs)
        if not np.any(finite):
            return 1

        best_comp = 1
        best_dist = np.inf
        for c in range(self.n_components):
            ref = self.means[c, vnir]
            diff = np.where(finite, obs - ref, 0.0)
            dist = np.sum(diff ** 2)
            if dist < best_dist:
                best_dist = dist
                best_comp = c

        return best_comp

    def classify_image(self, reflectance_stack):
        """Classify all pixels in an image stack.

        Args:
            reflectance_stack: 3D array [n_bands, rows, cols].

        Returns:
            2D array [rows, cols] of component indices.
        """
        n_bands, rows, cols = reflectance_stack.shape
        vnir = self.vnir_mask

        # Reshape to [n_bands, n_pixels]
        pixels = reflectance_stack.reshape(n_bands, -1)  # [n_bands, n_pix]
        vnir_pixels = pixels[vnir, :]  # [n_vnir, n_pix]

        # Compute distance to each component
        comp_labels = np.ones(pixels.shape[1], dtype=int)  # default soil
        best_dist = np.full(pixels.shape[1], np.inf)

        for c in range(self.n_components):
            ref = self.means[c, vnir][:, np.newaxis]  # [n_vnir, 1]
            diff = vnir_pixels - ref
            diff = np.where(np.isfinite(diff), diff, 0.0)
            dist = np.sum(diff ** 2, axis=0)  # [n_pix]
            better = dist < best_dist
            comp_labels[better] = c
            best_dist[better] = dist[better]

        return comp_labels.reshape(rows, cols)

    def prior_mean(self, component_idx):
        """Return prior mean spectrum for a component.

        Args:
            component_idx: Integer component index.

        Returns:
            1D array [n_bands] of mean reflectance.
        """
        return self.means[component_idx]

    def prior_variance(self, component_idx):
        """Return per-band variance for a component.

        Args:
            component_idx: Integer component index.

        Returns:
            1D array [n_bands] of variance values.
        """
        return self.variances[component_idx]

    def regularize(self, reflectance_stack, sigma_obs2=None, weight=0.1):
        """MAP regularization: blend observed reflectance with surface prior.

        For each pixel:
          1. Classify to nearest surface component
          2. Compute MAP estimate per band (diagonal covariance):
             r_map = (r_obs/sigma_obs^2 + r_prior/sigma_prior^2) /
                     (1/sigma_obs^2 + 1/sigma_prior^2)

        When sigma_obs2 is not provided, uses weight to control prior strength:
          sigma_obs^2 = prior_variance * (1 - weight) / weight

        Args:
            reflectance_stack: 3D array [n_bands, rows, cols].
            sigma_obs2: Optional 3D array [n_bands, rows, cols] of measurement
                        variance. If provided, uses actual uncertainty for MAP.
            weight: Prior weight when sigma_obs2 is not provided (0-1, default 0.1).

        Returns:
            3D array [n_bands, rows, cols] of regularized reflectance.
        """
        reflectance_stack = np.asarray(reflectance_stack, dtype=float)
        n_bands, rows, cols = reflectance_stack.shape
        result = reflectance_stack.copy()

        # Classify all pixels
        comp_labels = self.classify_image(reflectance_stack)

        for c in range(self.n_components):
            mask = comp_labels == c
            if not np.any(mask):
                continue

            xa = self.means[c]       # [n_bands]
            sa = self.variances[c]   # [n_bands]

            # Extract pixels for this component
            pix_data = reflectance_stack[:, mask]  # [n_bands, n_pix]

            if sigma_obs2 is not None:
                pix_sigma = sigma_obs2[:, mask]  # [n_bands, n_pix]
            else:
                # Use weight to set effective measurement variance
                pix_sigma = sa[:, np.newaxis] * (1.0 - weight) / max(weight, 1e-10)

            # MAP estimate (per-band, diagonal)
            xa_exp = xa[:, np.newaxis]   # [n_bands, 1]
            sa_exp = sa[:, np.newaxis]   # [n_bands, 1]

            inv_so = 1.0 / np.maximum(pix_sigma, 1e-20)
            inv_sa = 1.0 / np.maximum(sa_exp, 1e-20)

            pix_map = (pix_data * inv_so + xa_exp * inv_sa) / (inv_so + inv_sa)

            # Only regularize where observation is finite
            finite = np.isfinite(pix_data)
            pix_result = np.where(finite, pix_map, pix_data)

            result[:, mask] = pix_result

        return result
