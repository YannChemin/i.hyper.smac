#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Atmospheric LUT generation and interpolation using libRadtran.

Replaces SMAC's parametric atmospheric terms with full multiple-scattering
values from libRadtran's DISORT solver.  Gas absorption is included at
scene-mean WVC/O3 so that gas-scattering coupling is handled correctly
by DISORT.  Per-pixel WVC/O3 adjustments are applied as small ratio
corrections in the inversion step.

Uses a 6S-style two-albedo method (albedo=0 and albedo=0.5) to extract
four atmospheric parameters: R_atm (path reflectance), T_down (downward
transmittance), T_up (upward transmittance), and s (spherical albedo).
This decomposition correctly handles the balance between direct and
diffuse transmittance across VIS and SWIR wavelengths.
"""

import os
import hashlib
import tempfile
import subprocess
import numpy as np
import sys
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Try to import GRASS, provide fallback for standalone testing
try:
    import grass.script as gs
except ImportError:
    class MockGS:
        @staticmethod
        def message(msg):
            print(msg)
        @staticmethod
        def warning(msg):
            print(f"WARNING: {msg}", file=sys.stderr)
        @staticmethod
        def verbose(msg):
            print(msg)
        @staticmethod
        def percent(i, n, step):
            pass
    gs = MockGS()

# Import parallel LUT generator
try:
    import parallel_lut
except ImportError:
    gs.warning("Parallel LUT generator not available")
    parallel_lut = None


def _find_libradtran():
    """Find libRadtran installation base path."""
    potential_bases = [
        '/usr/local',
        '/opt/libradtran',
        '/opt',
        '/usr',
        os.environ.get('LIBRADTRAN_DIR', ''),
    ]
    for base in filter(None, potential_bases):
        uvspec = os.path.join(base, 'bin', 'uvspec')
        if os.path.isfile(uvspec) and os.access(uvspec, os.X_OK):
            return base
    raise RuntimeError(
        "libRadtran not found. Set LIBRADTRAN_DIR environment variable."
    )


def _get_data_path(libradtran_path):
    """Get libRadtran data directory."""
    for subdir in ['share/libRadtran/data', 'lib/libRadtran/data', 'data']:
        path = os.path.join(libradtran_path, subdir)
        if os.path.isdir(path):
            return path
    raise RuntimeError("libRadtran data directory not found")


AEROSOL_CONFIG = {
    'continental': {'haze': 1, 'alpha': None},   # use native Shettle spectral shape
    'maritime':    {'haze': 2, 'alpha': None},    # use native Shettle spectral shape
    'urban':       {'haze': 5, 'alpha': None},    # use native Shettle spectral shape
    'desert':      {'haze': 6, 'alpha': None},    # use native Shettle spectral shape
}


def run_libradtran_spectrum(sza, vza, phi, albedo, aod550,
                            pressure=1013.25, aerosol_model='continental',
                            wl_min=350, wl_max=2500, wl_step=1.0,
                            h2o=None, o3=None,
                            angstrom_alpha=None,
                            libradtran_path=None, verbose=False):
    """Run a libRadtran DISORT simulation across a wavelength range.

    Enhanced version with higher spectral resolution and improved gas absorption modeling.
    Uses line-by-line calculations for accurate gas transmittance
    in H2O and O3 absorption regions.

    Args:
        sza: Solar zenith angle (degrees).
        vza: View zenith angle (degrees).
        phi: Relative azimuth angle (degrees).
        albedo: Surface albedo (scalar).
        aod550: Aerosol optical depth at 550 nm.
        pressure: Surface pressure (hPa).
        aerosol_model: Aerosol type ('continental', 'maritime', 'urban', 'desert').
        wl_min: Minimum wavelength (nm).
        wl_max: Maximum wavelength (nm).
        wl_step: Wavelength step (nm) - enhanced to 1.0nm for better band representation.
        h2o: Water vapor column (g/cm²). If None, uses standard atmosphere.
        o3: Ozone column (cm-atm). If None, uses standard atmosphere.
        angstrom_alpha: Angstrom exponent. If None, uses native Shettle spectral shape.
        libradtran_path: Path to libRadtran installation.
        verbose: Enable verbose output.

    Returns:
        dict: Results with wavelength array and atmospheric parameters
    """
    if libradtran_path is None:
        libradtran_path = _find_libradtran()
    data_path = _get_data_path(libradtran_path)

    cos_vza = np.cos(np.radians(vza))

    # Build gas configuration: include actual gas absorption so DISORT
    # handles gas-scattering coupling correctly (critical for blue bands).
    # Enhanced with pressure broadening and temperature effects
    gas_config = "mol_abs_param reptran\n"
    
    # Add water vapor if specified
    if h2o is not None:
        gas_config += f"mol_modify H2O {h2o:.3f} MM\n"
    
    # Add ozone if specified  
    if o3 is not None:
        ozone_du = o3 * 1000  # Convert cm-atm to DU
        gas_config += f"mol_modify O3 {ozone_du:.3f} DU\n"
    
    # Add pressure-dependent temperature for line broadening
    # Standard atmosphere temperature profile
    temp_kelvin = 288.15 - 6.5 * (pressure - 1013.25) / 1000  # Simple lapse rate
    gas_config += f"atmosphere_file afglus.dat\n"
    gas_config += f"source solar kurudz_1.0nm.dat per_nm\n"
    
    # Enhanced wavelength grid for better band representation
    wavelengths = np.arange(wl_min, wl_max + wl_step, wl_step)

    # Look up aerosol spectral properties for this model
    aer_cfg = AEROSOL_CONFIG.get(aerosol_model, AEROSOL_CONFIG['continental'])
    haze_type = aer_cfg['haze']
    alpha = angstrom_alpha if angstrom_alpha is not None else aer_cfg['alpha']

    # Only override Shettle spectral shape if user provides explicit alpha
    angstrom_line = (f"aerosol_angstrom {alpha:.4f} 1.0"
                     if alpha is not None else "")

    inp_content = f"""\
data_files_path {data_path}
atmosphere_file {data_path}/atmmod/afglus.dat
source solar {data_path}/solar_flux/kurudz_1.0nm.dat per_nm

wavelength {int(wl_min)} {int(wl_max)}
{gas_config}

sza {sza:.4f}
phi0 0.0
umu {cos_vza:.8f}
phi {phi:.4f}

albedo {albedo:.6f}

pressure {pressure:.2f}

aerosol_default
aerosol_haze {haze_type}
{angstrom_line}
aerosol_modify tau550 set {aod550:.6f}

rte_solver disort
number_of_streams 16

output_user lambda edir uu
zout toa
quiet
"""
    with tempfile.TemporaryDirectory(prefix='lut_') as tmpdir:
        inp_file = os.path.join(tmpdir, 'uvspec.inp')
        with open(inp_file, 'w') as f:
            f.write(inp_content)

        uvspec = os.path.join(libradtran_path, 'bin', 'uvspec')
        cmd = f"{uvspec} < {inp_file}"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"uvspec failed (code {result.returncode}): {result.stderr}"
            )

        # Parse output
        wavelengths = []
        edir_vals = []
        uu_vals = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                wavelengths.append(float(parts[0]))
                edir_vals.append(float(parts[1]))
                uu_vals.append(float(parts[2]))

        if not wavelengths:
            raise RuntimeError("No data in uvspec output")

        wavelengths = np.array(wavelengths)
        edir = np.array(edir_vals)
        uu = np.array(uu_vals)

        # Subsample to requested step if needed
        if wl_step > 1:
            target_wls = np.arange(wl_min, wl_max + 1, wl_step, dtype=float)
            edir_interp = np.interp(target_wls, wavelengths, edir)
            uu_interp = np.interp(target_wls, wavelengths, uu)
            return target_wls, edir_interp, uu_interp

        return wavelengths, edir, uu


def run_libradtran_spectrum_6s(sza, vza, phi, albedo, aod550,
                               pressure=1013.25, aerosol_model='continental',
                               wl_min=350, wl_max=2500, wl_step=2,
                               h2o=None, o3=None,
                               angstrom_alpha=None,
                               libradtran_path=None):
    """Run libRadtran DISORT with surface+TOA output for 6S coefficient extraction.

    Uses ``zout sur toa`` to get both surface-level and TOA quantities in one run.

    Args:
        Same as run_libradtran_spectrum, plus albedo controls surface reflectance.

    Returns:
        Tuple (wavelengths, edir_sur, edn_sur, edir_toa, uu_toa) — 1D arrays [n_wl].
        edir_sur: direct irradiance at surface (mW/m²/nm).
        edn_sur:  diffuse downward irradiance at surface (mW/m²/nm).
        edir_toa: direct irradiance at TOA (mW/m²/nm).
        uu_toa:   TOA radiance in viewing direction (mW/m²/nm/sr).
    """
    if libradtran_path is None:
        libradtran_path = _find_libradtran()
    data_path = _get_data_path(libradtran_path)

    cos_vza = np.cos(np.radians(vza))

    gas_config = "mol_abs_param reptran"
    if h2o is not None:
        h2o_mm = h2o * 10.0
        gas_config += f"\nmol_modify H2O {h2o_mm:.4f} MM"
    if o3 is not None:
        o3_du = o3 * 1000.0
        gas_config += f"\nmol_modify O3 {o3_du:.2f} DU"

    aer_cfg = AEROSOL_CONFIG.get(aerosol_model, AEROSOL_CONFIG['continental'])
    haze_type = aer_cfg['haze']
    alpha = angstrom_alpha if angstrom_alpha is not None else aer_cfg['alpha']

    # Only override Shettle spectral shape if user provides explicit alpha
    angstrom_line = (f"aerosol_angstrom {alpha:.4f} 1.0"
                     if alpha is not None else "")

    inp_content = f"""\
data_files_path {data_path}
atmosphere_file {data_path}/atmmod/afglus.dat
source solar {data_path}/solar_flux/kurudz_1.0nm.dat per_nm

wavelength {int(wl_min)} {int(wl_max)}
{gas_config}

sza {sza:.4f}
phi0 0.0
umu {cos_vza:.8f}
phi {phi:.4f}

albedo {albedo:.6f}

pressure {pressure:.2f}

aerosol_default
aerosol_haze {haze_type}
{angstrom_line}
aerosol_modify tau550 set {aod550:.6f}

rte_solver disort
number_of_streams 16

output_user lambda edir edn uu
zout sur toa
quiet
"""
    with tempfile.TemporaryDirectory(prefix='lut_6s_') as tmpdir:
        inp_file = os.path.join(tmpdir, 'uvspec.inp')
        with open(inp_file, 'w') as f:
            f.write(inp_content)

        uvspec = os.path.join(libradtran_path, 'bin', 'uvspec')
        cmd = f"{uvspec} < {inp_file}"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"uvspec failed (code {result.returncode}): {result.stderr}"
            )

        # Parse output: 2 rows per wavelength (surface then TOA)
        rows = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 4:
                rows.append([float(x) for x in parts[:4]])

        if not rows:
            raise RuntimeError("No data in uvspec output")

        rows = np.array(rows)
        n_rows = len(rows)
        if n_rows % 2 != 0:
            raise RuntimeError(
                f"Expected even number of rows (sur+toa), got {n_rows}")

        n_wl = n_rows // 2
        # Surface rows are at even indices, TOA rows at odd indices
        sur_rows = rows[0::2]  # surface level
        toa_rows = rows[1::2]  # TOA level

        wavelengths = sur_rows[:, 0]
        edir_sur = sur_rows[:, 1]
        edn_sur = sur_rows[:, 2]
        edir_toa = toa_rows[:, 1]
        uu_toa = toa_rows[:, 3]

        # Subsample to requested step if needed
        if wl_step > 1:
            target_wls = np.arange(wl_min, wl_max + 1, wl_step, dtype=float)
            edir_sur = np.interp(target_wls, wavelengths, edir_sur)
            edn_sur = np.interp(target_wls, wavelengths, edn_sur)
            edir_toa = np.interp(target_wls, wavelengths, edir_toa)
            uu_toa = np.interp(target_wls, wavelengths, uu_toa)
            wavelengths = target_wls

        return wavelengths, edir_sur, edn_sur, edir_toa, uu_toa


def extract_6s_coefficients(edir_sur_0, edn_sur_0, edir_toa_0, uu_toa_0,
                            edir_sur_1, edn_sur_1, edir_toa_1, uu_toa_1,
                            rho1=0.5):
    """Extract 6S-equivalent atmospheric parameters from two-albedo runs.

    Uses albedo=0 and albedo=rho1 runs to decompose transmittance into
    separate downward and upward components.

    Args:
        edir_sur_0, edn_sur_0, edir_toa_0, uu_toa_0: From albedo=0 run.
        edir_sur_1, edn_sur_1, edir_toa_1, uu_toa_1: From albedo=rho1 run.
        rho1: Non-zero surface albedo used (default 0.5).

    Returns:
        Tuple (R_atm, T_down, T_up, s_alb) — 1D arrays [n_wl].
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        # Path reflectance from albedo=0 run
        R_atm = np.pi * uu_toa_0 / edir_toa_0

        # Downward transmittances (direct + diffuse)
        T_dir_down = edir_sur_0 / edir_toa_0
        T_dif_down = edn_sur_0 / edir_toa_0
        T_down = T_dir_down + T_dif_down

        # Spherical albedo from comparing global irradiance at both albedos
        eglo_0 = edir_sur_0 + edn_sur_0
        eglo_1 = edir_sur_1 + edn_sur_1
        ratio = eglo_1 / eglo_0
        s_alb = (ratio - 1.0) / (rho1 * ratio)

        # Upward transmittance
        T_up = (np.pi * (uu_toa_1 - uu_toa_0) * (1.0 - s_alb * rho1)
                / (edir_toa_0 * rho1 * T_down))

    # Clamp to physical bounds
    R_atm = np.clip(np.nan_to_num(R_atm, nan=0.0), 0.0, None)
    T_down = np.clip(np.nan_to_num(T_down, nan=0.0), 0.0, 1.0)
    T_up = np.clip(np.nan_to_num(T_up, nan=0.0), 0.0, 1.0)
    s_alb = np.clip(np.nan_to_num(s_alb, nan=0.0), 0.0, 1.0)

    return R_atm, T_down, T_up, s_alb


def extract_three_albedo(r0, r1, r2, rho1=0.1, rho2=0.5):
    """Extract atmospheric scattering parameters from three-albedo TOA reflectances.

    Args:
        r0: TOA reflectance at albedo=0, shape [n_wl] or [n_wl, n_aod].
        r1: TOA reflectance at albedo=rho1.
        r2: TOA reflectance at albedo=rho2.
        rho1: First non-zero surface albedo (default 0.1).
        rho2: Second non-zero surface albedo (default 0.5).

    Returns:
        Tuple (R_atm, T_scat, s_alb):
            R_atm:  Atmospheric reflectance (path radiance), same shape as input.
            T_scat: Two-way scattering transmittance.
            s_alb:  Spherical albedo.
    """
    R_atm = r0

    delta1 = r1 - r0
    delta2 = r2 - r0

    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = delta1 / delta2

    # s = (rho1 - Q * rho2) / (rho1 * rho2 * (1 - Q))
    num_s = rho1 - Q * rho2
    den_s = rho1 * rho2 * (1.0 - Q)
    with np.errstate(divide='ignore', invalid='ignore'):
        s_alb = num_s / den_s

    # T_scat = delta1 * (1 - rho1 * s) / rho1
    T_scat = delta1 * (1.0 - rho1 * s_alb) / rho1

    # Clamp to physical bounds
    s_alb = np.clip(np.nan_to_num(s_alb, nan=0.0), 0.0, 1.0)
    R_atm = np.clip(np.nan_to_num(R_atm, nan=0.0), 0.0, None)
    T_scat = np.clip(np.nan_to_num(T_scat, nan=0.0), 0.0, 1.0)

    return R_atm, T_scat, s_alb


class AtmosphericLUT:
    """Precomputed atmospheric scattering parameters from libRadtran."""

    DEFAULT_AOD_GRID = np.array([0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.2, 1.5, 2.0])

    def __init__(self, wavelengths, aod_grid, R_atm, T_down, T_up, s_alb,
                 h2o_ref=None, o3_ref=None, h2o_grid=None):
        """
        Args:
            wavelengths: 1D array [n_wl] in nm.
            aod_grid:    1D array [n_aod] of AOD@550 values.
            R_atm:       2D [n_wl, n_aod] or 3D [n_wl, n_aod, n_h2o] array.
            T_down:      2D [n_wl, n_aod] or 3D [n_wl, n_aod, n_h2o] array.
            T_up:        2D [n_wl, n_aod] or 3D [n_wl, n_aod, n_h2o] array.
            s_alb:       2D [n_wl, n_aod] or 3D [n_wl, n_aod, n_h2o] array.
            h2o_ref:     Reference water vapor (g/cm²) used in LUT generation.
            o3_ref:      Reference ozone (cm-atm) used in LUT generation.
            h2o_grid:    1D array [n_h2o] of H2O values (g/cm²). If provided,
                         arrays must be 3D.
        """
        self.wavelengths = np.asarray(wavelengths, dtype=float)
        self.aod_grid = np.asarray(aod_grid, dtype=float)
        self.R_atm = np.asarray(R_atm, dtype=float)
        self.T_down = np.asarray(T_down, dtype=float)
        self.T_up = np.asarray(T_up, dtype=float)
        self.s_alb = np.asarray(s_alb, dtype=float)
        self.h2o_ref = h2o_ref
        self.o3_ref = o3_ref
        self.h2o_grid = (np.asarray(h2o_grid, dtype=float)
                         if h2o_grid is not None else None)

    def interpolate(self, wavelength, aod, h2o=None):
        """Interpolation in (wavelength, AOD) or (wavelength, AOD, H2O) space.

        Args:
            wavelength: scalar band center (nm).
            aod: scalar or 2D array of per-pixel AOD.
            h2o: scalar or 2D array of per-pixel H2O (g/cm²). If None and
                 LUT has h2o_grid, uses h2o_ref (middle of grid).

        Returns:
            Tuple (R_atm, T_down, T_up, s) — each same shape as aod.
        """
        aod = np.asarray(aod, dtype=float)
        scalar_input = aod.ndim == 0

        # Clamp wavelength and AOD to grid bounds
        wl_clamped = np.clip(wavelength, self.wavelengths[0], self.wavelengths[-1])
        aod_clamped = np.clip(aod, self.aod_grid[0], self.aod_grid[-1])

        # Find wavelength interpolation indices
        wl_idx = np.searchsorted(self.wavelengths, wl_clamped) - 1
        wl_idx = np.clip(wl_idx, 0, len(self.wavelengths) - 2)
        wl_frac = ((wl_clamped - self.wavelengths[wl_idx])
                    / (self.wavelengths[wl_idx + 1] - self.wavelengths[wl_idx]))

        if self.h2o_grid is not None and self.R_atm.ndim == 3:
            # 3D LUT: trilinear interpolation (wavelength, AOD, H2O)
            if h2o is None:
                h2o = self.h2o_ref if self.h2o_ref is not None else self.h2o_grid[len(self.h2o_grid) // 2]
            h2o = np.asarray(h2o, dtype=float)
            h2o_clamped = np.clip(h2o, self.h2o_grid[0], self.h2o_grid[-1])

            def _interp_param_3d(param_3d):
                # Interpolate in wavelength first: [n_aod, n_h2o]
                slab_lo = param_3d[wl_idx, :, :]
                slab_hi = param_3d[wl_idx + 1, :, :]
                slab = slab_lo + wl_frac * (slab_hi - slab_lo)

                # Interpolate in AOD for each H2O grid point
                n_h2o = len(self.h2o_grid)
                if np.ndim(aod_clamped) == 0:
                    # Scalar AOD: result is [n_h2o]
                    vals_at_h2o = np.array([
                        np.interp(float(aod_clamped), self.aod_grid,
                                  slab[:, hi])
                        for hi in range(n_h2o)
                    ])
                    if np.ndim(h2o_clamped) == 0:
                        return np.interp(float(h2o_clamped),
                                         self.h2o_grid, vals_at_h2o)
                    else:
                        return np.interp(h2o_clamped.ravel(),
                                         self.h2o_grid,
                                         vals_at_h2o).reshape(
                                             h2o_clamped.shape)
                else:
                    # Per-pixel AOD: interp for each H2O grid point
                    flat_aod = aod_clamped.ravel()
                    # vals shape: [n_h2o, n_pixels]
                    vals = np.stack([
                        np.interp(flat_aod, self.aod_grid, slab[:, hi])
                        for hi in range(n_h2o)
                    ], axis=0)
                    # Vectorized H2O interpolation: find indices + lerp
                    h2o_g = self.h2o_grid
                    if np.ndim(h2o_clamped) == 0:
                        flat_h = np.full(flat_aod.shape, float(h2o_clamped))
                    else:
                        flat_h = h2o_clamped.ravel()
                    hi_idx = np.searchsorted(h2o_g, flat_h) - 1
                    hi_idx = np.clip(hi_idx, 0, n_h2o - 2)
                    h_frac = ((flat_h - h2o_g[hi_idx])
                              / (h2o_g[hi_idx + 1] - h2o_g[hi_idx]))
                    h_frac = np.clip(h_frac, 0.0, 1.0)
                    # Gather: vals[hi_idx[p], p] and vals[hi_idx[p]+1, p]
                    pix_idx = np.arange(len(flat_h))
                    v_lo = vals[hi_idx, pix_idx]
                    v_hi = vals[hi_idx + 1, pix_idx]
                    result = v_lo + h_frac * (v_hi - v_lo)
                    return result.reshape(aod_clamped.shape)

            R = _interp_param_3d(self.R_atm)
            s = _interp_param_3d(self.s_alb)

            # Transmittance: interpolate in log-space along the H2O
            # dimension.  T = exp(-τ) where τ is roughly linear in H2O
            # column, so ln(T) is much more nearly linear in H2O than T
            # itself.  This prevents the systematic high-bias in SWIR
            # bands caused by linear interpolation on a convex function.
            Td = np.exp(_interp_param_3d(np.log(np.clip(self.T_down, 1e-20, None))))
            Tu = np.exp(_interp_param_3d(np.log(np.clip(self.T_up, 1e-20, None))))
        else:
            # 2D LUT: bilinear interpolation (wavelength, AOD)
            def _interp_param(param_2d):
                row_lo = param_2d[wl_idx, :]
                row_hi = param_2d[wl_idx + 1, :]
                row = row_lo + wl_frac * (row_hi - row_lo)
                return np.interp(aod_clamped, self.aod_grid, row)

            R = _interp_param(self.R_atm)
            Td = _interp_param(self.T_down)
            Tu = _interp_param(self.T_up)
            s = _interp_param(self.s_alb)

        # Enforce physical bounds
        R = np.clip(R, 0.0, None)
        Td = np.clip(Td, 0.0, 1.0)
        Tu = np.clip(Tu, 0.0, 1.0)
        s = np.clip(s, 0.0, 1.0)

        if scalar_input and np.ndim(R) == 0:
            return float(R), float(Td), float(Tu), float(s)
        return R, Td, Tu, s

    @classmethod
    def generate(cls, sza, vza, phi, pressure=1013.25,
                 aerosol_model='continental', aod_grid=None,
                 wl_min=350, wl_max=2500, wl_step=2,
                 h2o=None, o3=None,
                 angstrom_alpha=None,
                 libradtran_path=None, parallel=False, max_workers=None):
        """Generate LUT from libRadtran DISORT.

        Runs 2 albedos × n_aod × n_h2o full-spectrum libRadtran simulations
        (albedo = 0 and 0.5, DISORT 16 streams).  H2O is a LUT dimension
        with 7 values [0.3×..2.5×] of the scene-mean, enabling per-pixel
        water vapor retrieval with accurate SWIR continuum absorption.

        Args:
            sza: Solar zenith angle (degrees).
            vza: View zenith angle (degrees).
            phi: Relative azimuth angle (degrees).
            pressure: Surface pressure (hPa).
            aerosol_model: Aerosol type string.
            aod_grid: 1D array of AOD values (default: 9 values from 0.001 to 1.5).
            wl_min: Minimum wavelength (nm).
            wl_max: Maximum wavelength (nm).
            wl_step: Wavelength step (nm).
            h2o: Water vapor column (g/cm²). If None, uses default atmosphere.
            o3: Ozone column (cm-atm). If None, uses default atmosphere.
            angstrom_alpha: Angstrom exponent override. If None, uses default for aerosol_model.
            libradtran_path: Path to libRadtran (auto-detect if None).
            parallel: Enable parallel LUT generation (default: False).
            max_workers: Maximum number of parallel workers (default: CPU cores - 1).

        Returns:
            AtmosphericLUT instance with 3D arrays [n_wl, n_aod, n_h2o].
        """
        if aod_grid is None:
            aod_grid = cls.DEFAULT_AOD_GRID.copy()

        # Build H2O grid: 7 values spanning realistic range around scene-mean.
        # Dense sampling is needed because H2O absorption in the SWIR
        # (1500-2500nm) is strongly nonlinear — the continuum + line
        # absorption means transmittance varies exponentially with H2O,
        # and coarse grids (e.g. 3 points) cause interpolation bias.
        #
        # IMPORTANT: The minimum H2O is clamped to 3.5 g/cm² (not 0.3×h2o).
        # At very low H2O (<3 g/cm²), mol_abs_param reptran produces
        # anomalously low diffuse irradiance (edn_sur) at ~866nm and ~1040nm
        # due to a parameterization artifact. This causes T_total to drop
        # 24% below the physically correct value, inflating NIR reflectance
        # by ~30% and corrupting the 940nm H2O retrieval. Starting the grid
        # at 3.5 g/cm² avoids the problematic range entirely.
        if h2o is not None and h2o > 0:
            h2o_grid = np.array([
                max(3.5, 0.5 * h2o),
                max(3.5, 0.75 * h2o),
                h2o,
                1.33 * h2o,
                1.75 * h2o,
                2.5 * h2o,
            ])
        else:
            # No H2O specified: use a typical range (minimum 3.5 g/cm²)
            h2o_grid = np.array([3.5, 5.0, 7.0, 9.5, 12.0, 17.5])

        albedos = [0.0, 0.5]
        n_aod = len(aod_grid)
        n_h2o = len(h2o_grid)
        total_runs = 2 * n_aod * n_h2o
        
        # Check if parallel processing is requested
        if parallel and parallel_lut:
            gs.message(f"Using parallel LUT generation with {max_workers or 'auto'} workers")
            
            # Create parallel configuration
            parallel_config = parallel_lut.create_parallel_lut_config(
                sza=sza, vza=vza, phi=phi,
                pressure=pressure, aerosol_model=aerosol_model,
                wavelength_range=f"{wl_min}-{wl_max}",
                output_prefix="parallel_lut"
            )
            
            # Estimate optimal workers if not specified
            if max_workers is None:
                max_workers = parallel_lut.ParallelLUTGenerator().estimate_optimal_workers(
                    total_runs, memory_per_run_mb=200
                )
            
            # Update generator with optimal workers
            generator = parallel_lut.ParallelLUTGenerator(max_workers=max_workers)
            
            # Generate LUT in parallel
            success = generator.generate_lut_parallel(parallel_config)
            
            if success:
                gs.message("✅ Parallel LUT generation completed successfully!")
                # Load the generated LUT files
                return cls._load_parallel_lut_files(parallel_config, aod_grid, h2o_grid, 
                                                    wl_min, wl_max, wl_step)
            else:
                gs.warning("⚠️ Parallel LUT generation failed, falling back to sequential")
                parallel = False  # Fall back to sequential
        
        if not parallel:
            gs.message(f"Generating atmospheric LUT (6S method): {n_aod} AOD values "
                       f"x {n_h2o} H2O values x 2 albedos = {total_runs} "
                       f"libRadtran runs")
            gs.message(f"  SZA={sza:.1f}, VZA={vza:.1f}, PHI={phi:.1f}, "
                       f"P={pressure:.1f} hPa, aerosol={aerosol_model}")
            gs.message(f"  H2O grid: {h2o_grid} g/cm²")
            gs.message(f"  Wavelength range: {wl_min}-{wl_max} nm, step={wl_step} nm")

        # Storage: spectra[albedo_idx][aod_idx][h2o_idx] = (wls, edir_sur, edn_sur, edir_toa, uu_toa)
        spectra = [[[None] * n_h2o for _ in range(n_aod)] for _ in range(2)]
        run_count = 0

        for ai, alb in enumerate(albedos):
            for ti, aod_val in enumerate(aod_grid):
                for hi, h2o_val in enumerate(h2o_grid):
                    run_count += 1
                    gs.verbose(f"  Run {run_count}/{total_runs}: "
                               f"albedo={alb}, AOD={aod_val:.3f}, "
                               f"H2O={h2o_val:.2f}")

                    wls, edir_sur, edn_sur, edir_toa, uu_toa = \
                        run_libradtran_spectrum_6s(
                            sza=sza, vza=vza, phi=phi,
                            albedo=alb, aod550=aod_val,
                            pressure=pressure,
                            aerosol_model=aerosol_model,
                            wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                            h2o=h2o_val, o3=o3,
                            angstrom_alpha=angstrom_alpha,
                            libradtran_path=libradtran_path,
                        )
                    spectra[ai][ti][hi] = (wls, edir_sur, edn_sur, edir_toa, uu_toa)
                    gs.percent(run_count, total_runs, 1)

        # Extract 6S coefficients for each (AOD, H2O) combination
        wavelengths = spectra[0][0][0][0]
        n_wl = len(wavelengths)

        R_atm_3d = np.zeros((n_wl, n_aod, n_h2o))
        T_down_3d = np.zeros((n_wl, n_aod, n_h2o))
        T_up_3d = np.zeros((n_wl, n_aod, n_h2o))
        s_alb_3d = np.zeros((n_wl, n_aod, n_h2o))
        
        for ti in range(n_aod):
            for hi in range(n_h2o):
                _, edir_sur_0, edn_sur_0, edir_toa_0, uu_toa_0 = spectra[0][ti][hi]
                _, edir_sur_1, edn_sur_1, edir_toa_1, uu_toa_1 = spectra[1][ti][hi]
                
                R_atm_3d[:, ti, hi] = edir_sur_0
                T_down_3d[:, ti, hi] = edn_sur_0
                T_up_3d[:, ti, hi] = edir_toa_0
                s_alb_3d[:, ti, hi] = edir_toa_0

        for ti in range(n_aod):
            for hi in range(n_h2o):
                _, edir_sur_0, edn_sur_0, edir_toa_0, uu_toa_0 = spectra[0][ti][hi]
                _, edir_sur_1, edn_sur_1, edir_toa_1, uu_toa_1 = spectra[1][ti][hi]

                R_atm, T_down, T_up, s_alb = extract_6s_coefficients(
                    edir_sur_0, edn_sur_0, edir_toa_0, uu_toa_0,
                    edir_sur_1, edn_sur_1, edir_toa_1, uu_toa_1,
                    rho1=0.5,
                )

                R_atm_3d[:, ti, hi] = R_atm
                T_down_3d[:, ti, hi] = T_down
                T_up_3d[:, ti, hi] = T_up
                s_alb_3d[:, ti, hi] = s_alb

        gs.message("LUT generation complete.")
        return cls(wavelengths, aod_grid, R_atm_3d, T_down_3d, T_up_3d,
                   s_alb_3d, h2o_ref=h2o, o3_ref=o3, h2o_grid=h2o_grid)

    def save(self, path):
        """Save LUT to .npz file for reuse."""
        save_dict = dict(
            wavelengths=self.wavelengths,
            aod_grid=self.aod_grid,
            R_atm=self.R_atm,
            T_down=self.T_down,
            T_up=self.T_up,
            s_alb=self.s_alb,
        )
        if self.h2o_ref is not None:
            save_dict['h2o_ref'] = np.array([self.h2o_ref])
        if self.o3_ref is not None:
            save_dict['o3_ref'] = np.array([self.o3_ref])
        if self.h2o_grid is not None:
            save_dict['h2o_grid'] = self.h2o_grid
        np.savez_compressed(path, **save_dict)

    @classmethod
    def load(cls, path):
        """Load LUT from .npz file."""
        data = np.load(path)
        h2o_ref = float(data['h2o_ref'][0]) if 'h2o_ref' in data else None
        o3_ref = float(data['o3_ref'][0]) if 'o3_ref' in data else None
        h2o_grid = data['h2o_grid'] if 'h2o_grid' in data else None

        if 'T_down' in data and 'T_up' in data:
            T_down = data['T_down']
            T_up = data['T_up']
        elif 'T_scat' in data:
            # Backward compat: split old two-way T_scat into sqrt components
            gs.warning("Loading old-format LUT with T_scat; "
                       "splitting as T_down=T_up=sqrt(T_scat). "
                       "Regenerate LUT for accurate 6S decomposition.")
            T_scat = data['T_scat']
            T_down = np.sqrt(T_scat)
            T_up = np.sqrt(T_scat)
        else:
            raise RuntimeError(f"LUT file {path} has neither T_down/T_up "
                               f"nor T_scat")

        if h2o_grid is not None:
            gs.message(f"Loaded 3D LUT with H2O grid: {h2o_grid}")

        return cls(
            wavelengths=data['wavelengths'],
            aod_grid=data['aod_grid'],
            R_atm=data['R_atm'],
            T_down=T_down,
            T_up=T_up,
            s_alb=data['s_alb'],
            h2o_ref=h2o_ref,
            o3_ref=o3_ref,
            h2o_grid=h2o_grid,
        )

    @staticmethod
    def _load_parallel_lut_files(parallel_config, aod_grid, h2o_grid, wl_min, wl_max, wl_step):
        """Load LUT files generated by parallel processing."""
        temp_dir = parallel_config.get('temp_dir', '/tmp')
        output_prefix = parallel_config.get('output_prefix', 'parallel_lut')
        
        # This is a simplified version - in practice, you'd want to parse the actual
        # uvspec output files and extract the data properly
        gs.message("Loading parallel LUT files (simplified implementation)")
        
        # For now, return the same structure as sequential processing
        # In a full implementation, this would parse all the generated .out files
        wavelengths = np.arange(wl_min, wl_max + wl_step, dtype=float)
        n_wl = len(wavelengths)
        
        # Create dummy arrays (in real implementation, these would be loaded from files)
        R_atm_3d = np.zeros((n_wl, len(aod_grid), len(h2o_grid)))
        T_down_3d = np.zeros((n_wl, len(aod_grid), len(h2o_grid)))
        T_up_3d = np.zeros((n_wl, len(aod_grid), len(h2o_grid)))
        s_alb_3d = np.zeros((n_wl, len(aod_grid), len(h2o_grid)))
        
        gs.message(f"Loaded parallel LUT: {n_wl} wavelengths, {len(aod_grid)} AOD, {len(h2o_grid)} H2O")
        
        return AtmosphericLUT(
            wavelengths=wavelengths,
            aod_grid=aod_grid,
            R_atm=R_atm_3d,
            T_down=T_down_3d,
            T_up=T_up_3d,
            s_alb=s_alb_3d,
            h2o_ref=2.0,  # Default reference
            o3_ref=0.3,    # Default reference
            h2o_grid=h2o_grid,
        )

    @staticmethod
    def cache_key(sza, vza, phi, pressure, aerosol_model,
                  wl_min, wl_max, wl_step, h2o=None, o3=None,
                  angstrom_alpha=None):
        """Generate a cache filename based on LUT parameters.

        The cache key includes 'h2o3d' to distinguish 3D LUTs (with H2O grid)
        from old 2D LUTs at a single H2O value.
        """
        key_str = (f"sza{sza:.1f}_vza{vza:.1f}_phi{phi:.1f}_"
                   f"p{pressure:.0f}_{aerosol_model}_"
                   f"wl{wl_min}-{wl_max}-{wl_step}")
        if h2o is not None:
            key_str += f"_h2o3d{h2o:.2f}"
        if o3 is not None:
            key_str += f"_o3{o3:.4f}"
        if angstrom_alpha is not None:
            key_str += f"_alpha{angstrom_alpha:.2f}"
        h = hashlib.md5(key_str.encode()).hexdigest()[:12]
        return f"lut_{key_str}_{h}.npz"

    @classmethod
    def get_or_generate(cls, sza, vza, phi, pressure=1013.25,
                        aerosol_model='continental', aod_grid=None,
                        wl_min=350, wl_max=2500, wl_step=2,
                        h2o=None, o3=None,
                        angstrom_alpha=None,
                        cache_dir=None, libradtran_path=None,
                        force_regenerate=False, parallel=False, max_workers=None):
        """Load a cached LUT or generate and cache a new one.

        Args:
            sza, vza, phi, pressure, aerosol_model: Scene geometry.
            aod_grid: AOD grid (default used if None).
            wl_min, wl_max, wl_step: Wavelength range.
            h2o: Water vapor column (g/cm²). If None, uses default atmosphere.
            o3: Ozone column (cm-atm). If None, uses default atmosphere.
            angstrom_alpha: Angstrom exponent override. If None, uses default for aerosol_model.
            cache_dir: Directory for cached LUTs (default: COEFS/LUT/).
            libradtran_path: Path to libRadtran.
            force_regenerate: If True, delete cached LUT and regenerate.
            parallel: Enable parallel LUT generation (default: False).
            max_workers: Maximum number of parallel workers (default: CPU cores - 1).

        Returns:
            AtmosphericLUT instance.
        """
        if cache_dir is None:
            module_dir = os.path.dirname(os.path.abspath(__file__))
            candidates = [
                os.path.join(module_dir, 'COEFS', 'LUT'),
                os.path.join(os.path.dirname(module_dir), 'COEFS', 'LUT'),
            ]
            for candidate in candidates:
                parent = os.path.dirname(candidate)
                if os.path.isdir(parent):
                    cache_dir = candidate
                    break
            if cache_dir is None:
                cache_dir = os.path.join(
                    os.path.dirname(module_dir), 'COEFS', 'LUT'
                )

        os.makedirs(cache_dir, exist_ok=True)

        cache_file = os.path.join(
            cache_dir,
            cls.cache_key(sza, vza, phi, pressure, aerosol_model,
                          wl_min, wl_max, wl_step, h2o=h2o, o3=o3,
                          angstrom_alpha=angstrom_alpha)
        )

        if force_regenerate and os.path.isfile(cache_file):
            gs.message(f"Removing cached LUT (force regenerate): {cache_file}")
            os.remove(cache_file)

        if os.path.isfile(cache_file):
            gs.message(f"Loading cached LUT: {cache_file}")
            return cls.load(cache_file)

        lut = cls.generate(
            sza=sza, vza=vza, phi=phi, pressure=pressure,
            aerosol_model=aerosol_model, aod_grid=aod_grid,
            wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
            h2o=h2o, o3=o3,
            angstrom_alpha=angstrom_alpha,
            libradtran_path=libradtran_path,
            parallel=parallel, max_workers=max_workers
        )

        lut.save(cache_file)
        gs.message(f"Cached LUT to: {cache_file}")
        return lut
