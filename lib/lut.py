#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Atmospheric LUT generation and interpolation using libRadtran.

Replaces SMAC's single-scattering parametric terms (R_atm, T_scat, s) with
full multiple-scattering values from libRadtran's DISORT solver.  Gas
absorption is handled separately via SMAC coefficients, so the LUT only
contains scattering parameters (Rayleigh + aerosol, gas-free).

The three-albedo method extracts R_atm, T_scat, s from three libRadtran
runs at surface albedos 0, 0.1, and 0.5.
"""

import os
import hashlib
import tempfile
import subprocess
import numpy as np
import sys

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


AEROSOL_TYPE_MAP = {
    'continental': 'continental_average',
    'maritime': 'maritime_clean',
    'urban': 'urban',
    'desert': 'desert',
}


def run_libradtran_spectrum(sza, vza, phi, albedo, aod550,
                            pressure=1013.25, aerosol_model='continental',
                            wl_min=350, wl_max=2500, wl_step=10,
                            libradtran_path=None):
    """Run a single gas-free libRadtran DISORT simulation across a wavelength range.

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
        wl_step: Wavelength step (nm).
        libradtran_path: Path to libRadtran installation (auto-detect if None).

    Returns:
        Tuple (wavelengths, edir, uu) — 1D arrays of shape [n_wl].
        edir is direct solar irradiance at TOA (mW/m^2/nm),
        uu is TOA radiance in viewing direction (mW/m^2/nm/sr).
    """
    if libradtran_path is None:
        libradtran_path = _find_libradtran()
    data_path = _get_data_path(libradtran_path)

    cos_vza = np.cos(np.radians(vza))

    # Build wavelength grid string for output
    # libRadtran with mol_abs_param none will output at each nm in the range
    inp_content = f"""\
data_files_path {data_path}
atmosphere_file {data_path}/atmmod/afglus.dat
source solar {data_path}/solar_flux/kurudz_1.0nm.dat per_nm

wavelength {int(wl_min)} {int(wl_max)}
mol_abs_param none

sza {sza:.4f}
phi0 0.0
umu {cos_vza:.8f}
phi {phi:.4f}

albedo {albedo:.6f}

pressure {pressure:.2f}

aerosol_default
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

    DEFAULT_AOD_GRID = np.array([0.001, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5])

    def __init__(self, wavelengths, aod_grid, R_atm, T_scat, s_alb):
        """
        Args:
            wavelengths: 1D array [n_wl] in nm.
            aod_grid:    1D array [n_aod] of AOD@550 values.
            R_atm:       2D array [n_wl, n_aod] atmospheric reflectance (gas-free).
            T_scat:      2D array [n_wl, n_aod] two-way scattering transmittance.
            s_alb:       2D array [n_wl, n_aod] spherical albedo.
        """
        self.wavelengths = np.asarray(wavelengths, dtype=float)
        self.aod_grid = np.asarray(aod_grid, dtype=float)
        self.R_atm = np.asarray(R_atm, dtype=float)
        self.T_scat = np.asarray(T_scat, dtype=float)
        self.s_alb = np.asarray(s_alb, dtype=float)

    def interpolate(self, wavelength, aod):
        """Bilinear interpolation in (wavelength, AOD) space.

        Args:
            wavelength: scalar band center (nm).
            aod: scalar or 2D array of per-pixel AOD.

        Returns:
            Tuple (R_atm, T_scat, s) — each same shape as aod.
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

        # For each parameter, interpolate in wavelength first, then AOD
        def _interp_param(param_2d):
            # Interpolate in wavelength: get two rows
            row_lo = param_2d[wl_idx, :]     # shape [n_aod]
            row_hi = param_2d[wl_idx + 1, :]  # shape [n_aod]
            row = row_lo + wl_frac * (row_hi - row_lo)  # shape [n_aod]

            # Now interpolate in AOD
            return np.interp(aod_clamped, self.aod_grid, row)

        R = _interp_param(self.R_atm)
        T = _interp_param(self.T_scat)
        s = _interp_param(self.s_alb)

        # Enforce physical bounds
        R = np.clip(R, 0.0, None)
        T = np.clip(T, 0.0, 1.0)
        s = np.clip(s, 0.0, 1.0)

        if scalar_input:
            return float(R), float(T), float(s)
        return R, T, s

    @classmethod
    def generate(cls, sza, vza, phi, pressure=1013.25,
                 aerosol_model='continental', aod_grid=None,
                 wl_min=350, wl_max=2500, wl_step=10,
                 libradtran_path=None):
        """Generate LUT from libRadtran DISORT.

        Runs 3 x len(aod_grid) full-spectrum libRadtran simulations
        (albedo = 0, 0.1, 0.5, gas-free, DISORT 16 streams).

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
            libradtran_path: Path to libRadtran (auto-detect if None).

        Returns:
            AtmosphericLUT instance.
        """
        if aod_grid is None:
            aod_grid = cls.DEFAULT_AOD_GRID.copy()

        albedos = [0.0, 0.1, 0.5]
        n_aod = len(aod_grid)
        total_runs = 3 * n_aod

        gs.message(f"Generating atmospheric LUT: {n_aod} AOD values x 3 albedos "
                   f"= {total_runs} libRadtran runs")
        gs.message(f"  SZA={sza:.1f}, VZA={vza:.1f}, PHI={phi:.1f}, "
                   f"P={pressure:.1f} hPa, aerosol={aerosol_model}")
        gs.message(f"  Wavelength range: {wl_min}-{wl_max} nm, step={wl_step} nm")

        # Storage: spectra[albedo_idx][aod_idx] = (wavelengths, edir, uu)
        spectra = [[None] * n_aod for _ in range(3)]
        run_count = 0

        for ai, alb in enumerate(albedos):
            for ti, aod_val in enumerate(aod_grid):
                run_count += 1
                gs.verbose(f"  Run {run_count}/{total_runs}: "
                           f"albedo={alb}, AOD={aod_val:.3f}")

                wls, edir, uu = run_libradtran_spectrum(
                    sza=sza, vza=vza, phi=phi,
                    albedo=alb, aod550=aod_val,
                    pressure=pressure,
                    aerosol_model=aerosol_model,
                    wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                    libradtran_path=libradtran_path,
                )
                spectra[ai][ti] = (wls, edir, uu)
                gs.percent(run_count, total_runs, 1)

        # Convert radiance to TOA reflectance: r = pi * uu / edir
        # and extract R_atm, T_scat, s for each (wavelength, aod) pair
        wavelengths = spectra[0][0][0]  # reference wavelength grid
        n_wl = len(wavelengths)

        R_atm_2d = np.zeros((n_wl, n_aod))
        T_scat_2d = np.zeros((n_wl, n_aod))
        s_alb_2d = np.zeros((n_wl, n_aod))

        for ti in range(n_aod):
            # Compute TOA reflectance for each albedo
            _, edir0, uu0 = spectra[0][ti]
            _, edir1, uu1 = spectra[1][ti]
            _, edir2, uu2 = spectra[2][ti]

            # Use edir from albedo=0 run as the reference solar irradiance
            # (edir at TOA is independent of surface albedo)
            with np.errstate(divide='ignore', invalid='ignore'):
                r0 = np.pi * uu0 / edir0
                r1 = np.pi * uu1 / edir1
                r2 = np.pi * uu2 / edir2

            r0 = np.nan_to_num(r0, nan=0.0)
            r1 = np.nan_to_num(r1, nan=0.0)
            r2 = np.nan_to_num(r2, nan=0.0)

            R_atm, T_scat, s_alb = extract_three_albedo(
                r0, r1, r2, rho1=0.1, rho2=0.5
            )

            R_atm_2d[:, ti] = R_atm
            T_scat_2d[:, ti] = T_scat
            s_alb_2d[:, ti] = s_alb

        gs.message("LUT generation complete.")
        return cls(wavelengths, aod_grid, R_atm_2d, T_scat_2d, s_alb_2d)

    def save(self, path):
        """Save LUT to .npz file for reuse."""
        np.savez_compressed(
            path,
            wavelengths=self.wavelengths,
            aod_grid=self.aod_grid,
            R_atm=self.R_atm,
            T_scat=self.T_scat,
            s_alb=self.s_alb,
        )

    @classmethod
    def load(cls, path):
        """Load LUT from .npz file."""
        data = np.load(path)
        return cls(
            wavelengths=data['wavelengths'],
            aod_grid=data['aod_grid'],
            R_atm=data['R_atm'],
            T_scat=data['T_scat'],
            s_alb=data['s_alb'],
        )

    @staticmethod
    def cache_key(sza, vza, phi, pressure, aerosol_model,
                  wl_min, wl_max, wl_step):
        """Generate a cache filename based on LUT parameters."""
        key_str = (f"sza{sza:.1f}_vza{vza:.1f}_phi{phi:.1f}_"
                   f"p{pressure:.0f}_{aerosol_model}_"
                   f"wl{wl_min}-{wl_max}-{wl_step}")
        h = hashlib.md5(key_str.encode()).hexdigest()[:12]
        return f"lut_{key_str}_{h}.npz"

    @classmethod
    def get_or_generate(cls, sza, vza, phi, pressure=1013.25,
                        aerosol_model='continental', aod_grid=None,
                        wl_min=350, wl_max=2500, wl_step=10,
                        cache_dir=None, libradtran_path=None):
        """Load a cached LUT or generate and cache a new one.

        Args:
            sza, vza, phi, pressure, aerosol_model: Scene geometry.
            aod_grid: AOD grid (default used if None).
            wl_min, wl_max, wl_step: Wavelength range.
            cache_dir: Directory for cached LUTs (default: COEFS/LUT/).
            libradtran_path: Path to libRadtran.

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
                          wl_min, wl_max, wl_step)
        )

        if os.path.isfile(cache_file):
            gs.message(f"Loading cached LUT: {cache_file}")
            return cls.load(cache_file)

        lut = cls.generate(
            sza=sza, vza=vza, phi=phi, pressure=pressure,
            aerosol_model=aerosol_model, aod_grid=aod_grid,
            wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
            libradtran_path=libradtran_path,
        )

        lut.save(cache_file)
        gs.message(f"Cached LUT to: {cache_file}")
        return lut
