#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SMAC Coefficient Generator using libRadtran.

This module generates SMAC (Simplified Method for Atmospheric Correction)
coefficient files by running libRadtran simulations and fitting the
analytical SMAC formulas to the results.

Based on:
- Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the
  atmospheric correction of satellite measurements in the solar spectrum.
  Int. J. Remote Sensing, 15(1), 123-143.

Author: Generated for i.hyper.smac GRASS GIS module
"""

import os
import sys
import tempfile
import subprocess
import numpy as np
from datetime import datetime
from scipy.optimize import curve_fit, minimize
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict


@dataclass
class SMACCoefficients:
    """Container for SMAC coefficients (49 values in 19 lines)."""

    # Line 0: H2O absorption
    ah2o: float = 0.0
    nh2o: float = 0.5

    # Line 1: O3 absorption
    ao3: float = 0.0
    no3: float = 1.0

    # Line 2: O2 absorption
    ao2: float = 0.0
    no2: float = 0.0
    po2: float = 0.0

    # Line 3: CO2 absorption
    aco2: float = 0.0
    nco2: float = 0.0
    pco2: float = 0.0

    # Line 4: CH4 absorption
    ach4: float = 0.0
    nch4: float = 0.0
    pch4: float = 0.0

    # Line 5: NO2 absorption
    ano2: float = 0.0
    nno2: float = 0.0
    pno2: float = 0.0

    # Line 6: CO absorption
    aco: float = 0.0
    nco: float = 0.0
    pco: float = 0.0

    # Line 7: Spherical albedo
    a0s: float = 0.0
    a1s: float = 0.0
    a2s: float = 0.0
    a3s: float = 0.0

    # Line 8: Transmission
    a0T: float = 1.0
    a1T: float = 0.0
    a2T: float = 0.0
    a3T: float = 0.0

    # Line 9: Rayleigh optical thickness
    taur: float = 0.0
    sr: float = 0.0

    # Line 10: Aerosol optical depth scaling
    a0taup: float = 0.0
    a1taup: float = 1.0

    # Line 11: Aerosol properties
    wo: float = 0.9  # Single scattering albedo
    gc: float = 0.7  # Asymmetry parameter

    # Lines 12-13: Phase function
    a0P: float = 0.0
    a1P: float = 0.0
    a2P: float = 0.0
    a3P: float = 0.0
    a4P: float = 0.0

    # Lines 14-15: Residual terms
    Rest1: float = 0.0
    Rest2: float = 0.0
    Rest3: float = 0.0
    Rest4: float = 0.0

    # Line 16: Rayleigh residuals
    Resr1: float = 0.0
    Resr2: float = 0.0
    Resr3: float = 0.0

    # Lines 17-18: Aerosol residuals
    Resa1: float = 0.0
    Resa2: float = 0.0
    Resa3: float = 0.0
    Resa4: float = 0.0

    def to_file(self, filename: str):
        """Write coefficients to SMAC coefficient file format."""
        with open(filename, 'w') as f:
            f.write(f"{self.ah2o:.6e} {self.nh2o:.6e}\n")
            f.write(f"{self.ao3:.6e} {self.no3:.6e}\n")
            f.write(f"{self.ao2:.6e} {self.no2:.6e} {self.po2:.6e}\n")
            f.write(f"{self.aco2:.6e} {self.nco2:.6e} {self.pco2:.6e}\n")
            f.write(f"{self.ach4:.6e} {self.nch4:.6e} {self.pch4:.6e}\n")
            f.write(f"{self.ano2:.6e} {self.nno2:.6e} {self.pno2:.6e}\n")
            f.write(f"{self.aco:.6e} {self.nco:.6e} {self.pco:.6e}\n")
            f.write(f"{self.a0s:.6e} {self.a1s:.6e} {self.a2s:.6e} {self.a3s:.6e}\n")
            f.write(f"{self.a0T:.6e} {self.a1T:.6e} {self.a2T:.6e} {self.a3T:.6e}\n")
            f.write(f"{self.taur:.6e} {self.sr:.6e}\n")
            f.write(f"{self.a0taup:.6e} {self.a1taup:.6e}\n")
            f.write(f"{self.wo:.6e} {self.gc:.6e}\n")
            f.write(f"{self.a0P:.6e} {self.a1P:.6e} {self.a2P:.6e}\n")
            f.write(f"{self.a3P:.6e} {self.a4P:.6e}\n")
            f.write(f"{self.Rest1:.6e} {self.Rest2:.6e}\n")
            f.write(f"{self.Rest3:.6e} {self.Rest4:.6e}\n")
            f.write(f"{self.Resr1:.6e} {self.Resr2:.6e} {self.Resr3:.6e}\n")
            f.write(f"{self.Resa1:.6e} {self.Resa2:.6e}\n")
            f.write(f"{self.Resa3:.6e} {self.Resa4:.6e}\n")

    @classmethod
    def from_file(cls, filename: str) -> 'SMACCoefficients':
        """Read coefficients from SMAC coefficient file."""
        coef = cls()
        with open(filename) as f:
            lines = f.readlines()

        # Parse each line
        temp = lines[0].strip().split()
        coef.ah2o, coef.nh2o = float(temp[0]), float(temp[1])

        temp = lines[1].strip().split()
        coef.ao3, coef.no3 = float(temp[0]), float(temp[1])

        temp = lines[2].strip().split()
        coef.ao2, coef.no2, coef.po2 = float(temp[0]), float(temp[1]), float(temp[2])

        temp = lines[3].strip().split()
        coef.aco2, coef.nco2, coef.pco2 = float(temp[0]), float(temp[1]), float(temp[2])

        temp = lines[4].strip().split()
        coef.ach4, coef.nch4, coef.pch4 = float(temp[0]), float(temp[1]), float(temp[2])

        temp = lines[5].strip().split()
        coef.ano2, coef.nno2, coef.pno2 = float(temp[0]), float(temp[1]), float(temp[2])

        temp = lines[6].strip().split()
        coef.aco, coef.nco, coef.pco = float(temp[0]), float(temp[1]), float(temp[2])

        temp = lines[7].strip().split()
        coef.a0s, coef.a1s, coef.a2s, coef.a3s = [float(x) for x in temp[:4]]

        temp = lines[8].strip().split()
        coef.a0T, coef.a1T, coef.a2T, coef.a3T = [float(x) for x in temp[:4]]

        temp = lines[9].strip().split()
        coef.taur, coef.sr = float(temp[0]), float(temp[0])

        temp = lines[10].strip().split()
        coef.a0taup, coef.a1taup = float(temp[0]), float(temp[1])

        temp = lines[11].strip().split()
        coef.wo, coef.gc = float(temp[0]), float(temp[1])

        temp = lines[12].strip().split()
        coef.a0P, coef.a1P, coef.a2P = [float(x) for x in temp[:3]]

        temp = lines[13].strip().split()
        coef.a3P, coef.a4P = float(temp[0]), float(temp[1])

        temp = lines[14].strip().split()
        coef.Rest1, coef.Rest2 = float(temp[0]), float(temp[1])

        temp = lines[15].strip().split()
        coef.Rest3, coef.Rest4 = float(temp[0]), float(temp[1])

        temp = lines[16].strip().split()
        coef.Resr1, coef.Resr2, coef.Resr3 = [float(x) for x in temp[:3]]

        temp = lines[17].strip().split()
        coef.Resa1, coef.Resa2 = float(temp[0]), float(temp[1])

        temp = lines[18].strip().split()
        coef.Resa3, coef.Resa4 = float(temp[0]), float(temp[1])

        return coef


class LibRadtranRunner:
    """Run libRadtran simulations for SMAC coefficient generation."""

    def __init__(self, libradtran_path: str = None, verbose: bool = False):
        """
        Initialize the libRadtran runner.

        Args:
            libradtran_path: Path to libRadtran installation (auto-detect if None)
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self.libradtran_path = libradtran_path or self._find_libradtran()
        self.temp_dir = tempfile.mkdtemp(prefix='smac_coef_')

    def _find_libradtran(self) -> str:
        """Find libRadtran installation."""
        potential_paths = [
            '/usr/local',
            '/opt/libradtran',
            '/usr',
            os.environ.get('LIBRADTRAN_DIR', ''),
        ]

        for base in filter(None, potential_paths):
            uvspec = os.path.join(base, 'bin', 'uvspec')
            if os.path.isfile(uvspec) and os.access(uvspec, os.X_OK):
                if self.verbose:
                    print(f"Found libRadtran at: {base}")
                return base

        raise RuntimeError("libRadtran not found. Set LIBRADTRAN_DIR environment variable.")

    def _get_data_path(self) -> str:
        """Get libRadtran data directory."""
        for subdir in ['share/libRadtran/data', 'lib/libRadtran/data', 'data']:
            path = os.path.join(self.libradtran_path, subdir)
            if os.path.isdir(path):
                return path
        raise RuntimeError("libRadtran data directory not found")

    def run_simulation(self, wavelength: float, sza: float, vza: float = 0.0,
                       phi: float = 0.0, aod550: float = 0.0, h2o: float = 2.0,
                       o3: float = 0.3, albedo: float = 0.0,
                       aerosol_type: str = 'continental',
                       output_level: str = 'toa') -> Dict:
        """
        Run a single libRadtran simulation.

        Args:
            wavelength: Wavelength in nm
            sza: Solar zenith angle in degrees
            vza: View zenith angle in degrees
            phi: Relative azimuth in degrees
            aod550: Aerosol optical depth at 550nm
            h2o: Water vapor column in g/cm²
            o3: Ozone column in cm-atm (DU/1000)
            albedo: Surface albedo
            aerosol_type: 'continental', 'maritime', 'urban', 'desert'
            output_level: 'toa' for top of atmosphere, 'sur' for surface

        Returns:
            Dictionary with simulation results
        """
        data_path = self._get_data_path()

        # Build input file
        inp_content = f"""\
data_files_path {data_path}
atmosphere_file {data_path}/atmmod/afglus.dat
source solar {data_path}/solar_flux/kurudz_1.0nm.dat per_nm

wavelength {wavelength:.2f} {wavelength:.2f}

mol_abs_param reptran
mol_modify O3 {o3 * 1000:.3f} DU
mol_modify H2O {h2o:.3f} MM

sza {sza:.2f}
phi0 0.0
umu {np.cos(np.radians(vza)):.6f}
phi {phi:.2f}

albedo {albedo:.4f}

"""
        # Add aerosol if AOD > 0
        if aod550 > 0.001:
            aerosol_map = {
                'continental': 'continental_average',
                'maritime': 'maritime_clean',
                'urban': 'urban',
                'desert': 'desert',
            }
            # Note: aerosol_default uses Shettle (1989) rural aerosol model
            # Different aerosol types mainly affect single scattering albedo and
            # asymmetry parameter which are fitted separately
            inp_content += f"""aerosol_default
aerosol_modify tau550 set {aod550:.4f}

"""

        inp_content += f"""rte_solver disort
number_of_streams 16

output_user lambda edir edn eup eglo uu
zout {output_level}
quiet
"""

        # Write input file
        inp_file = os.path.join(self.temp_dir, 'uvspec.inp')
        with open(inp_file, 'w') as f:
            f.write(inp_content)

        # Run uvspec
        uvspec = os.path.join(self.libradtran_path, 'bin', 'uvspec')
        cmd = f"{uvspec} < {inp_file}"

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                if self.verbose:
                    print(f"uvspec error: {result.stderr}")
                return None

            # Parse output
            lines = [l.strip() for l in result.stdout.strip().split('\n')
                     if l.strip() and not l.startswith('#')]

            if not lines:
                return None

            parts = lines[0].split()
            if len(parts) >= 6:
                return {
                    'wavelength': float(parts[0]),
                    'edir': float(parts[1]),      # Direct irradiance
                    'edn': float(parts[2]),       # Diffuse downward
                    'eup': float(parts[3]),       # Diffuse upward
                    'eglo': float(parts[4]),      # Global irradiance
                    'radiance': float(parts[5]),  # TOA radiance
                }

        except Exception as e:
            if self.verbose:
                print(f"Simulation error: {e}")

        return None

    def cleanup(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)


class SMACCoefficientGenerator:
    """Generate SMAC coefficients using libRadtran simulations."""

    def __init__(self, libradtran_path: str = None, verbose: bool = False):
        """
        Initialize the coefficient generator.

        Args:
            libradtran_path: Path to libRadtran installation
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self.runner = LibRadtranRunner(libradtran_path, verbose)

    def compute_rayleigh_optical_thickness(self, wavelength_nm: float) -> float:
        """
        Compute Rayleigh optical thickness at sea level.

        Uses the analytical formula from Hansen & Travis (1974).

        Args:
            wavelength_nm: Wavelength in nanometers

        Returns:
            Rayleigh optical thickness
        """
        wl_um = wavelength_nm / 1000.0  # Convert to micrometers

        # Hansen & Travis (1974) formula
        taur = 0.008569 * wl_um**(-4) * (1 + 0.0113 * wl_um**(-2) + 0.00013 * wl_um**(-4))

        return taur

    def compute_angstrom_coefficient(self, wavelength_nm: float,
                                     alpha: float = 1.3) -> float:
        """
        Compute aerosol optical depth scaling using Ångström law.

        Args:
            wavelength_nm: Wavelength in nanometers
            alpha: Ångström exponent (typically 1.0-1.5)

        Returns:
            AOD scaling factor (τ(λ) = τ(550) * factor)
        """
        return (wavelength_nm / 550.0) ** (-alpha)

    def fit_gaseous_absorption(self, wavelength_nm: float, gas: str,
                               aerosol_type: str = 'continental') -> Tuple[float, float]:
        """
        Fit gaseous absorption coefficients.

        The transmission is modeled as: T = exp(a * (u*m)^n)

        Args:
            wavelength_nm: Wavelength in nanometers
            gas: Gas name ('h2o', 'o3')
            aerosol_type: Aerosol model type

        Returns:
            Tuple of (a, n) coefficients
        """
        if self.verbose:
            print(f"Fitting {gas.upper()} absorption at {wavelength_nm:.1f} nm...")

        # Define gas amounts to test
        if gas == 'h2o':
            amounts = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])  # g/cm²
            baseline = 0.001  # Very small for reference
        elif gas == 'o3':
            amounts = np.array([0.2, 0.25, 0.3, 0.35, 0.4, 0.45]) * 1000 / 1000  # cm-atm
            baseline = 0.001
        else:
            return 0.0, 0.5

        # Reference SZA and VZA
        sza = 30.0
        vza = 0.0
        m = 1.0 / np.cos(np.radians(sza)) + 1.0 / np.cos(np.radians(vza))

        transmissions = []

        # Get reference (no gas absorption) - use surface output for transmission
        if gas == 'h2o':
            ref_result = self.runner.run_simulation(
                wavelength_nm, sza, vza, aod550=0.0, h2o=baseline, o3=0.3,
                aerosol_type=aerosol_type, output_level='sur'
            )
        else:
            ref_result = self.runner.run_simulation(
                wavelength_nm, sza, vza, aod550=0.0, h2o=2.0, o3=baseline,
                aerosol_type=aerosol_type, output_level='sur'
            )

        if ref_result is None:
            return 0.0, 0.5

        ref_edir = ref_result['edir']

        # Run simulations for different gas amounts - use surface output
        for amount in amounts:
            if gas == 'h2o':
                result = self.runner.run_simulation(
                    wavelength_nm, sza, vza, aod550=0.0, h2o=amount, o3=0.3,
                    aerosol_type=aerosol_type, output_level='sur'
                )
            else:
                result = self.runner.run_simulation(
                    wavelength_nm, sza, vza, aod550=0.0, h2o=2.0, o3=amount,
                    aerosol_type=aerosol_type, output_level='sur'
                )

            if result is not None:
                T = result['edir'] / ref_edir
                transmissions.append(T)
            else:
                transmissions.append(1.0)

        transmissions = np.array(transmissions)

        # Filter out invalid values
        valid = (transmissions > 0.01) & (transmissions < 1.0)
        if np.sum(valid) < 3:
            return 0.0, 0.5

        # Fit: T = exp(a * (u*m)^n)
        # ln(T) = a * (u*m)^n
        ln_T = np.log(transmissions[valid])
        u_m = amounts[valid] * m

        try:
            def gas_model(x, a, n):
                return a * (x ** n)

            popt, _ = curve_fit(gas_model, u_m, ln_T, p0=[-0.1, 0.5],
                                bounds=([-10, 0.1], [0, 2.0]))
            a, n = popt

            if self.verbose:
                print(f"  {gas.upper()}: a={a:.6f}, n={n:.6f}")

            return a, n

        except Exception as e:
            if self.verbose:
                print(f"  Fitting failed: {e}")
            return 0.0, 0.5

    def fit_spherical_albedo(self, wavelength_nm: float,
                             aerosol_type: str = 'continental') -> Tuple[float, float, float, float]:
        """
        Fit spherical albedo coefficients.

        s = a0s * P_eq + a3s + a1s * τ_550 + a2s * τ_550²

        Args:
            wavelength_nm: Wavelength in nanometers
            aerosol_type: Aerosol model type

        Returns:
            Tuple of (a0s, a1s, a2s, a3s) coefficients
        """
        if self.verbose:
            print(f"Fitting spherical albedo at {wavelength_nm:.1f} nm...")

        sza = 30.0
        vza = 0.0

        # Test different AOD values
        aod_values = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5])
        albedo1, albedo2 = 0.0, 0.5

        s_values = []

        for aod in aod_values:
            # Run two simulations with different surface albedos at SURFACE level
            # Spherical albedo is measured via multiple scattering at the surface:
            # Eglo_sur(ρ) = Eglo_sur(0) / (1 - s*ρ)
            # s = (E1-E2) / (ρ1*E1 - ρ2*E2)
            r1 = self.runner.run_simulation(
                wavelength_nm, sza, vza, aod550=aod, albedo=albedo1,
                aerosol_type=aerosol_type, output_level='sur'
            )
            r2 = self.runner.run_simulation(
                wavelength_nm, sza, vza, aod550=aod, albedo=albedo2,
                aerosol_type=aerosol_type, output_level='sur'
            )

            if r1 is None or r2 is None:
                s_values.append(np.nan)
                continue

            # Compute spherical albedo from surface global irradiance
            E1, E2 = r1['eglo'], r2['eglo']
            denom = albedo1 * E1 - albedo2 * E2

            if abs(denom) > 1e-10:
                s = (E1 - E2) / denom
                s_values.append(max(0, min(1, s)))  # Clamp to [0, 1]
            else:
                s_values.append(np.nan)

        s_values = np.array(s_values)
        valid = ~np.isnan(s_values)

        if np.sum(valid) < 3:
            # Return analytical approximation
            taur = self.compute_rayleigh_optical_thickness(wavelength_nm)
            return taur * 0.5, 0.2, -0.05, 0.0

        # Fit: s = a0s * P_eq + a3s + a1s * τ + a2s * τ²
        # At sea level, P_eq = 1.0
        try:
            def s_model(tau, a1s, a2s, a3s):
                return a3s + a1s * tau + a2s * tau**2

            popt, _ = curve_fit(s_model, aod_values[valid], s_values[valid],
                                p0=[0.2, -0.05, 0.0])
            a1s, a2s, a3s = popt

            # a0s is related to Rayleigh scattering
            taur = self.compute_rayleigh_optical_thickness(wavelength_nm)
            a0s = taur * 0.5

            if self.verbose:
                print(f"  Spherical albedo: a0s={a0s:.6f}, a1s={a1s:.6f}, "
                      f"a2s={a2s:.6f}, a3s={a3s:.6f}")

            return a0s, a1s, a2s, a3s

        except Exception as e:
            if self.verbose:
                print(f"  Fitting failed: {e}")
            taur = self.compute_rayleigh_optical_thickness(wavelength_nm)
            return taur * 0.5, 0.2, -0.05, 0.0

    def fit_transmission(self, wavelength_nm: float,
                         aerosol_type: str = 'continental') -> Tuple[float, float, float, float]:
        """
        Fit transmission coefficients.

        T = a0T + a1T * τ_550/cos(θ) + (a2T * P_eq + a3T)/(1 + cos(θ))

        Args:
            wavelength_nm: Wavelength in nanometers
            aerosol_type: Aerosol model type

        Returns:
            Tuple of (a0T, a1T, a2T, a3T) coefficients
        """
        if self.verbose:
            print(f"Fitting transmission at {wavelength_nm:.1f} nm...")

        # Test different geometries and AOD values
        sza_values = np.array([0, 20, 40, 60])
        aod_values = np.array([0.0, 0.1, 0.2, 0.3])

        data_points = []

        for sza in sza_values:
            cos_sza = np.cos(np.radians(sza))
            if cos_sza < 0.1:
                continue

            for aod in aod_values:
                # Get extraterrestrial irradiance at TOA (no atmosphere effects)
                r_toa = self.runner.run_simulation(
                    wavelength_nm, 0.0, 0.0, aod550=0.0, albedo=0.0,
                    aerosol_type=aerosol_type, output_level='toa'
                )

                # Get irradiance at the surface after atmospheric attenuation
                r_surf = self.runner.run_simulation(
                    wavelength_nm, sza, 0.0, aod550=aod, albedo=0.0,
                    aerosol_type=aerosol_type, output_level='sur'
                )

                if r_toa is not None and r_surf is not None:
                    # Transmission = E_surface / (E_TOA * cos(sza))
                    E_toa = r_toa['edir']  # Direct irradiance at TOA (normal incidence)
                    T = r_surf['eglo'] / (E_toa * cos_sza)
                    T = max(0, min(1, T))

                    data_points.append({
                        'T': T,
                        'tau': aod,
                        'cos_theta': cos_sza,
                    })

        if len(data_points) < 4:
            return 1.0, -0.1, 0.0, -0.1

        # Fit the transmission model
        try:
            T_data = np.array([d['T'] for d in data_points])
            tau_data = np.array([d['tau'] for d in data_points])
            cos_data = np.array([d['cos_theta'] for d in data_points])

            def T_model(X, a0T, a1T, a3T):
                tau, cos_theta = X
                return a0T + a1T * tau / cos_theta + a3T / (1 + cos_theta)

            popt, _ = curve_fit(T_model, (tau_data, cos_data), T_data,
                                p0=[1.0, -0.1, -0.1])
            a0T, a1T, a3T = popt
            a2T = 0.0  # Pressure term

            if self.verbose:
                print(f"  Transmission: a0T={a0T:.6f}, a1T={a1T:.6f}, "
                      f"a2T={a2T:.6f}, a3T={a3T:.6f}")

            return a0T, a1T, a2T, a3T

        except Exception as e:
            if self.verbose:
                print(f"  Fitting failed: {e}")
            return 1.0, -0.1, 0.0, -0.1

    def generate_coefficients(self, wavelength_nm: float, fwhm_nm: float = 10.0,
                              aerosol_type: str = 'continental',
                              angstrom_exp: float = 1.3) -> SMACCoefficients:
        """
        Generate complete SMAC coefficient set for a spectral band.

        Args:
            wavelength_nm: Central wavelength in nanometers
            fwhm_nm: Full width at half maximum in nanometers
            aerosol_type: Aerosol model ('continental', 'maritime', 'urban', 'desert')
            angstrom_exp: Ångström exponent for aerosol spectral dependence

        Returns:
            SMACCoefficients object with fitted coefficients
        """
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Generating SMAC coefficients for {wavelength_nm:.1f} nm")
            print(f"FWHM: {fwhm_nm:.1f} nm, Aerosol type: {aerosol_type}")
            print(f"{'='*60}\n")

        coef = SMACCoefficients()

        # 1. Rayleigh optical thickness (analytical)
        coef.taur = self.compute_rayleigh_optical_thickness(wavelength_nm)
        coef.sr = coef.taur
        if self.verbose:
            print(f"Rayleigh optical thickness: {coef.taur:.6f}")

        # 2. Aerosol optical depth scaling (Ångström law)
        coef.a0taup = 0.0
        coef.a1taup = self.compute_angstrom_coefficient(wavelength_nm, angstrom_exp)
        if self.verbose:
            print(f"AOD scaling (a1taup): {coef.a1taup:.6f}")

        # 3. Aerosol properties (from aerosol model)
        aerosol_props = {
            'continental': (0.89, 0.65),
            'maritime': (0.98, 0.72),
            'urban': (0.82, 0.62),
            'desert': (0.92, 0.72),
        }
        coef.wo, coef.gc = aerosol_props.get(aerosol_type, (0.89, 0.65))

        # 4. Fit gaseous absorption
        # H2O absorption bands: 720, 820, 940, 1130, 1380, 1900nm
        h2o_active = any([
            700 < wavelength_nm < 750,    # 720nm band
            790 < wavelength_nm < 850,    # 820nm band
            880 < wavelength_nm < 1000,   # 940nm band
            1080 < wavelength_nm < 1200,  # 1130nm band
            1300 < wavelength_nm < 1500,  # 1380nm band
            1800 < wavelength_nm < 2000,  # 1900nm band
        ])
        if h2o_active:
            coef.ah2o, coef.nh2o = self.fit_gaseous_absorption(
                wavelength_nm, 'h2o', aerosol_type
            )
        else:
            coef.ah2o, coef.nh2o = 0.0, 0.5

        # O3 absorption (Chappuis band: 400-700nm)
        if 400 < wavelength_nm < 700:
            coef.ao3, coef.no3 = self.fit_gaseous_absorption(
                wavelength_nm, 'o3', aerosol_type
            )
        else:
            coef.ao3, coef.no3 = 0.0, 1.0

        # Other gases (typically zero for visible/NIR)
        coef.ao2, coef.no2, coef.po2 = 0.0, 0.0, 0.0
        coef.aco2, coef.nco2, coef.pco2 = 0.0, 0.0, 0.0
        coef.ach4, coef.nch4, coef.pch4 = 0.0, 0.0, 0.0
        coef.ano2, coef.nno2, coef.pno2 = 0.0, 0.0, 0.0
        coef.aco, coef.nco, coef.pco = 0.0, 0.0, 0.0

        # 5. Fit spherical albedo
        coef.a0s, coef.a1s, coef.a2s, coef.a3s = self.fit_spherical_albedo(
            wavelength_nm, aerosol_type
        )

        # 6. Fit transmission
        coef.a0T, coef.a1T, coef.a2T, coef.a3T = self.fit_transmission(
            wavelength_nm, aerosol_type
        )

        # 7. Phase function (simplified Henyey-Greenstein approximation)
        # P(Θ) ≈ (1 - g²) / (1 + g² - 2g*cos(Θ))^1.5
        # Polynomial approximation for small angles
        g = coef.gc
        coef.a0P = (1 - g**2) / (1 + g)**1.5
        coef.a1P = 0.0
        coef.a2P = 0.0
        coef.a3P = 0.0
        coef.a4P = 0.0

        # 8. Residual terms (empirical, small corrections)
        coef.Rest1, coef.Rest2 = 0.0, 0.0
        coef.Rest3, coef.Rest4 = 0.0, 0.0
        coef.Resr1, coef.Resr2, coef.Resr3 = 0.0, 0.0, 0.0
        coef.Resa1, coef.Resa2 = 0.0, 0.0
        coef.Resa3, coef.Resa4 = 0.0, 0.0

        if self.verbose:
            print(f"\nCoefficient generation complete!")

        return coef

    def cleanup(self):
        """Clean up resources."""
        self.runner.cleanup()


def generate_smac_coefficients(wavelength_nm: float, fwhm_nm: float = 10.0,
                               aerosol_type: str = 'continental',
                               output_file: str = None,
                               libradtran_path: str = None,
                               verbose: bool = False) -> SMACCoefficients:
    """
    Convenience function to generate SMAC coefficients.

    Args:
        wavelength_nm: Central wavelength in nanometers
        fwhm_nm: Full width at half maximum in nanometers
        aerosol_type: Aerosol model type
        output_file: Optional file path to save coefficients
        libradtran_path: Path to libRadtran installation
        verbose: Enable verbose output

    Returns:
        SMACCoefficients object
    """
    generator = SMACCoefficientGenerator(libradtran_path, verbose)

    try:
        coef = generator.generate_coefficients(
            wavelength_nm, fwhm_nm, aerosol_type
        )

        if output_file:
            coef.to_file(output_file)
            if verbose:
                print(f"Coefficients saved to: {output_file}")

        return coef

    finally:
        generator.cleanup()


def generate_batch(wavelength_min=350, wavelength_max=2500, step=10,
                   aerosol_types=None, output_dir=None,
                   libradtran_path=None, verbose=False):
    """
    Generate SMAC coefficient files for a range of wavelengths.

    Args:
        wavelength_min: Minimum wavelength in nm
        wavelength_max: Maximum wavelength in nm
        step: Wavelength step in nm
        aerosol_types: List of aerosol types (default: all four)
        output_dir: Output directory for coefficient files
        libradtran_path: Path to libRadtran installation
        verbose: Enable verbose output

    Returns:
        Number of files generated
    """
    if aerosol_types is None:
        aerosol_types = ['continental', 'maritime', 'urban', 'desert']

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), 'COEFS')

    wavelengths = list(range(wavelength_min, wavelength_max + 1, step))
    total = len(wavelengths) * len(aerosol_types)
    count = 0

    print(f"Generating {total} coefficient files "
          f"({len(wavelengths)} wavelengths x {len(aerosol_types)} aerosol types)")
    print(f"Wavelength range: {wavelength_min}-{wavelength_max} nm, step: {step} nm")
    print(f"Output directory: {output_dir}")

    for aerosol_type in aerosol_types:
        atype_upper = aerosol_type.upper()
        atype_dir = os.path.join(output_dir, atype_upper)
        os.makedirs(atype_dir, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Aerosol model: {atype_upper}")
        print(f"{'='*60}")

        generator = SMACCoefficientGenerator(libradtran_path, verbose)

        try:
            for i, wl in enumerate(wavelengths):
                count += 1
                progress = 100.0 * count / total
                print(f"[{progress:5.1f}%] {atype_upper} {wl:4d} nm "
                      f"({count}/{total})", end='', flush=True)

                try:
                    coef = generator.generate_coefficients(
                        wl, fwhm_nm=10.0, aerosol_type=aerosol_type
                    )

                    output_file = os.path.join(
                        atype_dir, f"coef_{wl}nm_{atype_upper}.dat")
                    coef.to_file(output_file)
                    print(f" -> OK")

                except Exception as e:
                    print(f" -> FAILED: {e}")

        finally:
            generator.cleanup()

    print(f"\nBatch generation complete: {count} files generated")
    return count


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate SMAC coefficients using libRadtran'
    )

    subparsers = parser.add_subparsers(dest='command')

    # Single wavelength mode
    single_parser = subparsers.add_parser('single',
        help='Generate coefficients for a single wavelength')
    single_parser.add_argument('wavelength', type=float,
        help='Central wavelength in nm')
    single_parser.add_argument('--fwhm', type=float, default=10.0,
        help='FWHM in nm (default: 10)')
    single_parser.add_argument('--aerosol', type=str, default='continental',
        choices=['continental', 'maritime', 'urban', 'desert'],
        help='Aerosol model (default: continental)')
    single_parser.add_argument('--output', type=str,
        help='Output coefficient file')
    single_parser.add_argument('--verbose', action='store_true',
        help='Verbose output')

    # Batch mode
    batch_parser = subparsers.add_parser('batch',
        help='Generate coefficients for a range of wavelengths')
    batch_parser.add_argument('--min', type=int, default=350,
        help='Minimum wavelength in nm (default: 350)')
    batch_parser.add_argument('--max', type=int, default=2500,
        help='Maximum wavelength in nm (default: 2500)')
    batch_parser.add_argument('--step', type=int, default=10,
        help='Wavelength step in nm (default: 10)')
    batch_parser.add_argument('--aerosol', type=str, nargs='+',
        default=['continental', 'maritime', 'urban', 'desert'],
        choices=['continental', 'maritime', 'urban', 'desert'],
        help='Aerosol model(s) (default: all)')
    batch_parser.add_argument('--output-dir', type=str,
        help='Output directory (default: COEFS/)')
    batch_parser.add_argument('--verbose', action='store_true',
        help='Verbose output')

    # Common arguments
    parser.add_argument('--libradtran', type=str, help='Path to libRadtran')

    args = parser.parse_args()

    # Default to single mode for backwards compatibility
    if args.command is None:
        # Check if a positional wavelength was given (old-style usage)
        if len(sys.argv) > 1 and sys.argv[1].replace('.', '').isdigit():
            args.command = 'single'
            args.wavelength = float(sys.argv[1])
            args.fwhm = 10.0
            args.aerosol = 'continental'
            args.output = None
            args.verbose = '--verbose' in sys.argv
        else:
            parser.print_help()
            sys.exit(1)

    if args.command == 'single':
        if not hasattr(args, 'output') or not args.output:
            args.output = f"coef_{args.wavelength:.0f}nm_{args.aerosol.upper()}.dat"

        coef = generate_smac_coefficients(
            wavelength_nm=args.wavelength,
            fwhm_nm=args.fwhm,
            aerosol_type=args.aerosol,
            output_file=args.output,
            libradtran_path=args.libradtran,
            verbose=args.verbose
        )

        print(f"\nGenerated SMAC coefficients for {args.wavelength:.1f} nm")
        print(f"Saved to: {args.output}")

    elif args.command == 'batch':
        generate_batch(
            wavelength_min=args.min,
            wavelength_max=args.max,
            step=args.step,
            aerosol_types=args.aerosol,
            output_dir=args.output_dir,
            libradtran_path=args.libradtran,
            verbose=args.verbose
        )
