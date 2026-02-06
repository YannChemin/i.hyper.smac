#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
libRadtran interface for atmospheric correction parameter generation.

This module provides functionality to generate atmospheric correction parameters
using libRadtran for SMAC atmospheric correction.
"""

import os
import re
import tempfile
import subprocess
import numpy as np
import sys
from datetime import datetime

# Try to import GRASS, provide fallback for standalone testing
try:
    import grass.script as gs
except ImportError:
    # Fallback for standalone testing without GRASS
    class MockGS:
        @staticmethod
        def message(msg):
            print(msg)
        @staticmethod
        def warning(msg):
            print(f"WARNING: {msg}", file=sys.stderr)
        @staticmethod
        def error(msg):
            print(f"ERROR: {msg}", file=sys.stderr)
        @staticmethod
        def fatal(msg):
            print(f"FATAL: {msg}", file=sys.stderr)
            sys.exit(1)
    gs = MockGS()

# Import smac module for coefficient class
try:
    import smac
except ImportError:
    try:
        from . import smac
    except ImportError:
        gs.warning("Could not import smac module, some functions may not work")

class LibRadtranRunner:
    def __init__(self, verbose=False):
        """
        Initialize the LibRadtran runner.
        
        Args:
            verbose (bool): Enable verbose output
        """
        self.verbose = verbose
        self.temp_dir = tempfile.mkdtemp(prefix='lradtran_')
        self.libradtran_path = self._find_libradtran()
        
    def _find_libradtran(self):
        """Find the libRadtran installation path considering typical layouts."""
        # Common base paths
        potential_bases = [
            '/usr/local',
            '/opt',
            '/usr',
            os.environ.get('LIBRADTRAN_DIR', '')
        ]
        
        for base in filter(None, potential_bases):
            uvspec_path = os.path.join(base, 'bin', 'uvspec')
            lib_path = os.path.join(base, 'lib', 'libRadtran')
            share_path = os.path.join(base, 'share', 'libRadtran')
    
            # Check if uvspec executable exists
            if os.path.isfile(uvspec_path) and os.access(uvspec_path, os.X_OK):
                # Verify that at least one of the data/library directories exist
                if os.path.isdir(lib_path) or os.path.isdir(share_path):
                    if self.verbose:
                        gs.message(f"Found libRadtran at: {base}")
                    return base
    
        raise RuntimeError("libRadtran not found. Ensure uvspec is in PATH and data directories exist.")
    
    def _create_input_file(self, params, input_file, output_type='reflectance'):
        """
        Create a libRadtran input file.
        
        Args:
            params (dict): Dictionary containing the parameters for the simulation
            input_file (str): Path to the input file to create
            output_type (str): Type of output - 'reflectance', 'transmittance', or 'radiance'
        """
        with open(input_file, 'w') as f:
            dp = os.path.join(self.libradtran_path, 'share/libRadtran/data/')
            dp_atmosph = os.path.join(self.libradtran_path, 'share/libRadtran/data/atmmod/')
            
            f.write("rte_solver disort\n")
            f.write("number_of_streams 8\n")
            
            # Atmospheric parameters
            f.write(f"data_files_path {dp}\n")
            f.write(f"atmosphere_file {os.path.join(dp_atmosph, 'afglus.dat')}\n")
            
            # Solar and viewing geometry
            f.write(f"albedo {params.get('surface_albedo', 0.1)}\n")
            f.write(f"sza {params.get('solar_zenith', 30.0)}\n")
            f.write(f"phi0 {params.get('solar_azimuth', 0.0)}\n")
            
            # Pressure
            if 'pressure' in params:
                f.write(f"pressure {params['pressure']}\n")
            
            # Aerosols
            if 'aerosol_model' in params:
                f.write("aerosol_default\n")
                f.write(f"aerosol_modify tau550 set {params.get('aod_550', 0.1)}\n")
            
            # Molecular absorption
            f.write("mol_abs_param reptran\n")
            ozone_du = params.get('ozone', 0.3) * 1000  # Convert cm-atm to DU
            f.write(f"mol_modify O3 {ozone_du:.3f} DU\n")
            f.write(f"mol_modify H2O {params.get('water_vapor', 2.0):.3f} MM\n")
            
            # Wavelength
            f.write(f"wavelength {params['wavelength']}\n")

            # Output configuration based on type
            if output_type == 'transmittance':
                f.write("output_quantity transmittance\n")
            elif output_type == 'radiance':
                # For path radiance calculation
                f.write("zout 100.0\n")  # TOA altitude in km
                f.write("umu 1.0\n")     # Looking straight down
                f.write(f"phi {params.get('view_azimuth', 0.0)}\n")
                f.write("output_user lambda spher_alb uu\n")
            else:  # reflectance
                f.write("output_user lambda edir edn eup\n")
            
            f.write("quiet\n")
    
    def run_simulation(self, params, output_type='reflectance'):
        """
        Run libRadtran simulation.
        
        Args:
            params (dict): Dictionary containing the parameters for the simulation
            output_type (str): Type of output - 'reflectance', 'transmittance', or 'radiance'
            
        Returns:
            dict: Dictionary containing the simulation results
        """
        # Create input file
        input_file = os.path.join(self.temp_dir, f'lradtran_{output_type}.inp')
        self._create_input_file(params, input_file, output_type)
        
        # Set up command
        uvspec_path = os.path.join(self.libradtran_path, 'bin', 'uvspec')
        cmd = f"{uvspec_path} < {input_file}"
        
        if self.verbose:
            gs.message(f"Running: {cmd}")
        
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            
            return self._parse_output(result.stdout, output_type)
            
        except subprocess.CalledProcessError as e:
            gs.error(f"Error running libRadtran: {e}")
            gs.error(f"Stderr: {e.stderr}")
            raise
    
    def _parse_output(self, output, output_type):
        """
        Parse libRadtran output.
        
        Args:
            output (str): libRadtran output as string
            output_type (str): Type of output
            
        Returns:
            dict: Dictionary containing the parsed results
        """
        lines = output.strip().split('\n')
        if not lines or len(lines) < 1:
            raise ValueError("Invalid libRadtran output format")
            
        # Get last non-empty line (skip header)
        data_line = None
        for line in reversed(lines):
            if line.strip() and not line.startswith('#'):
                data_line = line
                break
        
        if not data_line:
            raise ValueError("No data found in libRadtran output")
            
        data = data_line.split()
            
        if output_type == 'transmittance':
            # Default output: wavelength edir edn eup ...
            return {
                'wavelength': float(data[0]),
                'direct_transmittance': float(data[1]),
                'diffuse_transmittance': float(data[2]),
                'total_transmittance': float(data[1]) + float(data[2]),
            }
        elif output_type == 'radiance':
            # output_user: lambda spher_alb uu
            return {
                'wavelength': float(data[0]),
                'spherical_albedo': float(data[1]),
                'path_radiance': float(data[2]),
            }
        else:  # reflectance
            # output_user: lambda edir edn eup
            return {
                'wavelength': float(data[0]),
                'direct_irradiance': float(data[1]),
                'diffuse_down': float(data[2]),
                'diffuse_up': float(data[3]),
            }
    
    def cleanup(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def __del__(self):
        """Destructor to ensure cleanup."""
        try:
            self.cleanup()
        except:
            pass

def earth_sun_distance(year, month, day):
    """
    Simplified Earth-Sun distance using day-of-year approximation.
    Accurate to ~0.5% (good enough for most reflectance calculations).
    """
    # Day of year (1-365/366)
    doy = datetime(year, month, day).timetuple().tm_yday
    
    # Earth orbital eccentricity approximation
    # Distance = 1 + 0.0167 * cos(2*pi*(doy-3)/365)
    beta = 2 * np.pi * (doy - 3) / 365.25
    dist = 1 + 0.01670963 * np.cos(beta) - 0.0000146 * np.cos(2 * beta)
    
    return dist  # in AU

def gaussian_rsp(wl, wl_center, fwhm):
    # Gaussian response normalized to max = 1
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-0.5 * ((wl - wl_center) / sigma)**2)

def E0(wl_center, fwhm,
       uvspec_bin="/usr/local/bin/uvspec",
       solar_file="/usr/local/share/libRadtran/data/solar_flux/kurudz_1.0nm.dat",
       atmosphere_file="/usr/local/share/libRadtran/data/atmmod/afglus.dat",
       verbose=True):
    """
    Compute band-integrated exo-atmospheric irradiance E0_band using libRadtran.

    Args:
        wl_center (float): Central wavelength in nm
        fwhm (float): Full width at half maximum in nm
        uvspec_bin (str): Path to uvspec executable
        solar_file (str): Path to solar flux data file
        atmosphere_file (str): Path to atmosphere model file
        verbose (bool): Enable verbose output

    Returns:
        float: Band-integrated exo-atmospheric irradiance, or None on error
    """
    if verbose:
        gs.message(f"\n{'='*80}")
        gs.message(f"Processing band: {wl_center:.2f} nm, FWHM: {fwhm:.2f} nm")
        gs.message(f"Using solar file: {solar_file}")
        gs.message(f"Using atmosphere file: {atmosphere_file}")
        gs.message(f"{'='*80}")

    # Calculate wavelength range (nm), clamp to integer boundaries
    # for compatibility with solar flux file grid
    wl_min = int(np.floor(wl_center - 2 * fwhm))
    wl_max = int(np.ceil(wl_center + 2 * fwhm))
    # Ensure within Kurucz solar file range (250-10000nm)
    wl_min = max(wl_min, 250)
    wl_max = min(wl_max, 10000)

    # Derive data_files_path from solar_file path
    data_path = os.path.dirname(os.path.dirname(solar_file))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Build uvspec input file for TOA extraterrestrial irradiance
        uvspec_inp = f"""\
data_files_path {data_path}
atmosphere_file {atmosphere_file}
source solar {solar_file} per_nm
wavelength {wl_min} {wl_max}
zout toa
output_user lambda edir
quiet
"""
        # Write input file
        inp_path = os.path.join(tmpdir, "uvspec.inp")
        with open(inp_path, 'w') as f:
            f.write(uvspec_inp)

        if verbose:
            gs.message(f"\nUVSPEC Input ({inp_path}):")
            gs.message("-" * 40)
            for line in uvspec_inp.split('\n'):
                if line.strip():
                    gs.message(line)
            gs.message("-" * 40)
            gs.message(f"Running uvspec for {wl_center:.2f} nm...")

        # Run uvspec with proper shell command (uvspec reads from stdin)
        cmd = f"{uvspec_bin} < {inp_path}"

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                shell=True
            )

            if result.stderr and verbose:
                gs.warning("\nUVSPEC Warnings/Errors:")
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        gs.warning(f"  {line}")

            if result.returncode != 0:
                gs.error(f"uvspec failed with code {result.returncode}")
                return None

            if not result.stdout.strip():
                gs.error("Empty output from uvspec")
                return None

        except Exception as e:
            gs.error(f"Error running uvspec: {e}")
            return None

        # Parse output directly from stdout
        try:
            lines = [l.strip() for l in result.stdout.strip().split('\n')
                     if l.strip() and not l.startswith('#')]

            if not lines:
                gs.error("No data lines in uvspec output")
                return None

            # Parse wavelength and irradiance data
            data = []
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        wl = float(parts[0])
                        irr = float(parts[1])
                        data.append([wl, irr])
                    except ValueError:
                        continue

            if not data:
                gs.error("Could not parse uvspec output data")
                return None

            data = np.array(data)
            lam_uv = data[:, 0]    # nm
            edir = data[:, 1]      # direct solar irradiance at TOA

            if verbose:
                gs.message(f"\nSuccessfully processed {wl_center:.2f} nm")
                gs.message(f"Output shape: {data.shape}")
                gs.message("First few data points:")
                for row in data[:3]:
                    gs.message(f"  {row[0]:.2f} nm: {row[1]:.3e}")

        except Exception as e:
            gs.error(f"Error processing uvspec output: {e}")
            return None

    # Interpolate to regular wavelength grid
    wavelengths = np.linspace(lam_uv.min(), lam_uv.max(), 1000)
    E0_lambda = np.interp(wavelengths, lam_uv, edir)

    # Build band response and integrate
    R = gaussian_rsp(wavelengths, wl_center, fwhm)
    num = np.trapezoid(E0_lambda * R, wavelengths)
    den = np.trapezoid(R, wavelengths)
    E0_band = num / den

    return E0_band

def find_coef_file(wavelength, aerosol_model='continental', coef_dir=None):
    """
    Find the nearest pre-generated SMAC coefficient file for a given wavelength.

    Args:
        wavelength (float): Central wavelength in nm
        aerosol_model (str): Aerosol model type (continental, maritime, urban, desert)
        coef_dir (str): Optional path to coefficient directory. If None, uses
                        the COEFS directory relative to this module.

    Returns:
        tuple: (coef_file_path, actual_wavelength) or (None, None) if not found
    """
    # Determine coefficient directory
    if coef_dir is None:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        # Search candidate locations:
        # 1. COEFS/ next to this module (installed in GRASS etc/i_hyper_lib/COEFS/)
        # 2. COEFS/ in the parent directory (development layout)
        candidates = [
            os.path.join(module_dir, 'COEFS'),
            os.path.join(os.path.dirname(module_dir), 'COEFS'),
        ]
        coef_dir = None
        for candidate in candidates:
            if os.path.isdir(candidate):
                coef_dir = candidate
                break
        if coef_dir is None:
            return None, None

    # Map aerosol model to directory name
    aerosol_dir_map = {
        'continental': 'CONTINENTAL',
        'maritime': 'MARITIME',
        'urban': 'URBAN',
        'desert': 'DESERT'
    }

    aerosol_dir_name = aerosol_dir_map.get(aerosol_model.lower(), 'CONTINENTAL')
    aerosol_dir = os.path.join(coef_dir, aerosol_dir_name)

    if not os.path.isdir(aerosol_dir):
        return None, None

    # Scan directory for available coefficient files
    pattern = re.compile(r'coef_(\d+)nm_' + aerosol_dir_name + r'\.dat')
    available_wavelengths = []
    for f in os.listdir(aerosol_dir):
        m = pattern.match(f)
        if m:
            available_wavelengths.append(int(m.group(1)))

    if not available_wavelengths:
        return None, None

    # Find nearest wavelength
    nearest_wl = min(available_wavelengths, key=lambda x: abs(x - wavelength))

    # Construct filename
    coef_filename = f"coef_{nearest_wl}nm_{aerosol_dir_name}.dat"
    coef_path = os.path.join(aerosol_dir, coef_filename)

    if os.path.isfile(coef_path):
        return coef_path, nearest_wl

    return None, None


def get_smac_parameters(wavelength, fwhm=10.0, sza=30.0, vza=0.0,
                       aod_550=0.1, water_vapor=2.0, ozone=0.3,
                       pressure=1013.25, surface_albedo=0.1,
                       aerosol_model='continental', coef_dir=None,
                       generate=False, verbose=False):
    """
    Get SMAC atmospheric correction parameters for a specific wavelength.

    By default, loads pre-generated SMAC coefficients from the COEFS directory.
    When generate=True, uses libRadtran and scipy to fit SMAC coefficients
    at runtime, saving them for future reuse.

    Args:
        wavelength (float): Central wavelength in nm
        fwhm (float): Full width at half maximum in nm
        sza (float): Solar zenith angle in degrees
        vza (float): View zenith angle in degrees
        aod_550 (float): Aerosol optical depth at 550nm
        water_vapor (float): Total column water vapor in g/cm²
        ozone (float): Total column ozone in cm-atm
        pressure (float): Atmospheric pressure in hPa
        surface_albedo (float): Surface albedo (0-1)
        aerosol_model (str): Aerosol model (continental, maritime, urban, desert)
        coef_dir (str): Optional path to coefficient directory
        generate (bool): If True, generate coefficients from libRadtran instead
                        of loading pre-generated files (requires libRadtran + scipy)
        verbose (bool): Enable verbose output

    Returns:
        smac.coeff: SMAC coefficient object for use with smac_inv

    Raises:
        FileNotFoundError: If no coefficient file is found and generate=False
        RuntimeError: If generation fails (libRadtran not found, scipy missing, etc.)
    """
    if generate:
        return _generate_and_cache_coefficients(
            wavelength, fwhm, aerosol_model, coef_dir, verbose
        )

    # Find the nearest coefficient file
    coef_path, nearest_wl = find_coef_file(wavelength, aerosol_model, coef_dir)

    if coef_path is None:
        raise FileNotFoundError(
            f"No SMAC coefficient file found for wavelength {wavelength} nm "
            f"and aerosol model '{aerosol_model}'. "
            f"Check that the COEFS directory exists and contains the required files. "
            f"Or use the -g flag to generate coefficients from libRadtran."
        )

    if verbose:
        wl_diff = abs(wavelength - nearest_wl)
        if wl_diff > 0:
            gs.message(f"Loading SMAC coefficients for {wavelength:.1f} nm "
                      f"(using nearest: {nearest_wl} nm, diff: {wl_diff:.1f} nm)")
        else:
            gs.message(f"Loading SMAC coefficients for {wavelength:.1f} nm")

    # Load and return the coefficients
    coefs = smac.coeff(coef_path)

    return coefs


def _generate_and_cache_coefficients(wavelength, fwhm, aerosol_model,
                                      coef_dir, verbose):
    """
    Generate SMAC coefficients from libRadtran and cache them to disk.

    Args:
        wavelength (float): Central wavelength in nm
        fwhm (float): Full width at half maximum in nm
        aerosol_model (str): Aerosol model type
        coef_dir (str): Optional coefficient directory for caching
        verbose (bool): Enable verbose output

    Returns:
        smac.coeff: SMAC coefficient object
    """
    try:
        from smac_coef_generator import generate_smac_coefficients
    except ImportError:
        try:
            from lib.smac_coef_generator import generate_smac_coefficients
        except ImportError:
            raise RuntimeError(
                "Cannot import smac_coef_generator. Make sure scipy is installed "
                "and smac_coef_generator.py is available."
            )

    # Round wavelength to nearest integer for file naming
    wl_int = int(round(wavelength))
    aerosol_upper = aerosol_model.upper()

    # Determine cache directory
    if coef_dir is None:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(module_dir, 'COEFS'),
            os.path.join(os.path.dirname(module_dir), 'COEFS'),
        ]
        for candidate in candidates:
            if os.path.isdir(candidate):
                coef_dir = candidate
                break
        if coef_dir is None:
            # Create COEFS next to this module
            coef_dir = os.path.join(os.path.dirname(module_dir), 'COEFS')

    # Check if already cached
    cache_dir = os.path.join(coef_dir, aerosol_upper)
    cache_file = os.path.join(cache_dir, f"coef_{wl_int}nm_{aerosol_upper}.dat")

    if os.path.isfile(cache_file):
        if verbose:
            gs.message(f"Using cached generated coefficients: {cache_file}")
        return smac.coeff(cache_file)

    # Generate coefficients
    if verbose:
        gs.message(f"Generating SMAC coefficients for {wavelength:.1f} nm "
                  f"({aerosol_model}) from libRadtran...")

    coef_obj = generate_smac_coefficients(
        wavelength_nm=wavelength,
        fwhm_nm=fwhm,
        aerosol_type=aerosol_model,
        verbose=verbose
    )

    # Save to cache
    os.makedirs(cache_dir, exist_ok=True)
    coef_obj.to_file(cache_file)
    if verbose:
        gs.message(f"Cached coefficients to: {cache_file}")

    # Load through smac.coeff for consistent return type
    return smac.coeff(cache_file)


def generate_smac_coefficients_from_libradtran(wavelength, fwhm=10.0,
                                               solar_zenith=30.0,
                                               solar_azimuth=0.0,
                                               view_zenith=0.0,
                                               view_azimuth=0.0,
                                               aod_550=0.1,
                                               water_vapor=2.0,
                                               ozone=0.3,
                                               pressure=1013.25,
                                               aerosol_model='continental',
                                               verbose=False):
    """
    EXPERIMENTAL: Generate SMAC coefficients from libRadtran simulations.

    WARNING: This function is not fully implemented and produces incorrect
    coefficients. Use get_smac_parameters() with pre-generated coefficient
    files instead.

    The proper implementation would require:
    1. Running multiple libRadtran simulations varying AOD, SZA, surface albedo
    2. Fitting the SMAC analytical formulas to the simulation results
    3. Extracting polynomial coefficients from the fits

    This is kept for future development but should NOT be used for production.

    Args:
        wavelength (float): Central wavelength in nm
        fwhm (float): Full width at half maximum in nm
        solar_zenith (float): Solar zenith angle in degrees
        solar_azimuth (float): Solar azimuth angle in degrees
        view_zenith (float): View zenith angle in degrees
        view_azimuth (float): View azimuth angle in degrees
        aod_550 (float): Aerosol optical depth at 550nm
        water_vapor (float): Water vapor content in g/cm² or mm
        ozone (float): Ozone content in cm-atm (will be converted to DU)
        pressure (float): Atmospheric pressure in hPa
        aerosol_model (str): Aerosol model type
        verbose (bool): Enable verbose output

    Returns:
        smac.coeff: SMAC coefficient object (WARNING: coefficients are approximate)
    """
    gs.warning("generate_smac_coefficients_from_libradtran() is experimental and "
               "produces approximate coefficients. Use get_smac_parameters() with "
               "pre-generated coefficient files for accurate results.")

    # Try to use pre-generated coefficients first
    coef_path, nearest_wl = find_coef_file(wavelength, aerosol_model)
    if coef_path is not None:
        if verbose:
            gs.message(f"Using pre-generated coefficients for {nearest_wl} nm")
        return smac.coeff(coef_path)

    # Fall back to analytical approximation (for wavelengths outside 400-2500nm range)
    if verbose:
        gs.message(f"Generating analytical SMAC coefficients for {wavelength} nm...")

    # Create temporary file for coefficients
    temp_dir = tempfile.mkdtemp(prefix='smac_coef_')
    temp_coef_file = os.path.join(temp_dir, f'smac_coef_{wavelength}.dat')

    try:
        with open(temp_coef_file, 'w') as f:
            # H2O absorption coefficients
            # Strong absorption bands: 720, 820, 940, 1130, 1380, 1900nm
            ah2o = 0.0
            if 700 < wavelength < 750 or 800 < wavelength < 850:
                ah2o = -0.02
            elif 900 < wavelength < 980:
                ah2o = -0.08
            elif 1100 < wavelength < 1160:
                ah2o = -0.05
            elif 1350 < wavelength < 1420:
                ah2o = -0.15
            elif 1850 < wavelength < 1950:
                ah2o = -0.20
            f.write(f"{ah2o:.6e} 5.000000e-01\n")

            # O3 absorption (Chappuis band 500-700nm)
            ao3 = 0.0
            if 500 < wavelength < 700:
                # Peak absorption around 600nm
                ao3 = -0.05 * np.exp(-((wavelength - 600) / 100)**2)
            f.write(f"{ao3:.6e} 1.000000e+00\n")

            # O2, CO2, CH4, NO2, CO - typically zero for most bands
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # O2
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # CO2
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # CH4
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # NO2
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # CO

            # Spherical albedo coefficients: s = a0s*P + a3s + a1s*tau + a2s*tau^2
            # Rayleigh contribution dominates at short wavelengths
            taur = 0.008569 * (wavelength / 1000) ** (-4) * (1 + 0.0113 * (wavelength/1000)**(-2))
            a0s = taur * 0.5  # Approximate spherical albedo from Rayleigh
            a1s = 0.2   # Aerosol contribution
            a2s = -0.05  # Second order
            a3s = 0.0
            f.write(f"{a0s:.6e} {a1s:.6e} {a2s:.6e} {a3s:.6e}\n")

            # Transmission coefficients: T = a0T + a1T*tau/cos(θ) + (a2T*P + a3T)/(1+cos(θ))
            # These are the critical coefficients for proper AOD-dependent correction
            a0T = 1.0     # Base transmission (clear atmosphere)
            a1T = -0.15   # AOD dependence (negative: more AOD = less transmission)
            a2T = 0.0     # Pressure dependence
            a3T = -0.10   # Geometric factor
            f.write(f"{a0T:.6e} {a1T:.6e} {a2T:.6e} {a3T:.6e}\n")

            # Rayleigh optical thickness
            f.write(f"{taur:.6e} {taur:.6e}\n")

            # Aerosol optical depth scaling: taup = a0taup + a1taup * taup550
            # Angstrom law: tau(λ) = tau(550) * (λ/550)^(-α)
            alpha = 1.3  # Typical Angstrom exponent for continental aerosol
            a0taup = 0.0
            a1taup = (wavelength / 550.0) ** (-alpha)
            f.write(f"{a0taup:.6e} {a1taup:.6e}\n")

            # Aerosol single scattering albedo and asymmetry parameter
            wo = 0.89  # Continental aerosol
            gc = 0.65  # Continental aerosol
            f.write(f"{wo:.6e} {gc:.6e}\n")

            # Phase function coefficients (Henyey-Greenstein approximation)
            # P(θ) = a0P + a1P*θ + a2P*θ² + a3P*θ³ + a4P*θ⁴
            a0P = (1 - gc**2) / (1 + gc)**1.5
            f.write(f"{a0P:.6e} 0.000000e+00 0.000000e+00\n")
            f.write("0.000000e+00 0.000000e+00\n")

            # Residual terms (small corrections)
            f.write("0.000000e+00 0.000000e+00\n")  # Rest1, Rest2
            f.write("0.000000e+00 0.000000e+00\n")  # Rest3, Rest4
            f.write("0.000000e+00 0.000000e+00 0.000000e+00\n")  # Resr1, Resr2, Resr3
            f.write("0.000000e+00 0.000000e+00\n")  # Resa1, Resa2
            f.write("0.000000e+00 0.000000e+00\n")  # Resa3, Resa4

        # Load the coefficients
        coefs = smac.coeff(temp_coef_file)
        return coefs

    finally:
        # Cleanup
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    # Example usage - can run standalone for testing
    gs.message("Testing within GRASS environment" if 'grass' in sys.modules else
               "Testing in standalone mode")

    print("=" * 60)
    print("Testing SMAC coefficient loading from pre-generated files")
    print("=" * 60)

    # Test multiple wavelengths
    test_wavelengths = [400, 550, 850, 1000, 1600, 2200]

    for wavelength in test_wavelengths:
        try:
            coefs = get_smac_parameters(
                wavelength=wavelength,
                aerosol_model='continental',
                verbose=True
            )

            print(f"\n--- Coefficients for {wavelength} nm ---")
            print(f"H2O: ah2o={coefs.ah2o:.6e}, nh2o={coefs.nh2o:.6f}")
            print(f"O3:  ao3={coefs.ao3:.6e}, no3={coefs.no3:.6f}")
            print(f"Rayleigh: taur={coefs.taur:.6f}")
            print(f"Aerosol scaling: a0taup={coefs.a0taup:.6f}, a1taup={coefs.a1taup:.6f}")
            print(f"Aerosol props: wo={coefs.wo:.6f}, gc={coefs.gc:.6f}")
            print(f"Transmission: a0T={coefs.a0T:.6f}, a1T={coefs.a1T:.6f}, a3T={coefs.a3T:.6f}")
            print(f"Sph. albedo: a0s={coefs.a0s:.6e}, a1s={coefs.a1s:.6f}")

        except Exception as e:
            print(f"Error loading coefficients for {wavelength} nm: {e}")

    # Test SMAC correction with proper coefficients
    print("\n" + "=" * 60)
    print("Testing SMAC atmospheric correction")
    print("=" * 60)

    try:
        coefs = get_smac_parameters(wavelength=850, aerosol_model='continental')

        # Test parameters
        r_toa = 0.27      # Typical NIR TOA reflectance
        theta_s = 30.0    # Solar zenith
        phi_s = 0.0
        theta_v = 0.0     # Nadir view
        phi_v = 0.0
        pressure = 1013.25
        uo3 = 0.3
        uh2o = 2.0

        print(f"\nInput: TOA reflectance = {r_toa}")
        print(f"       SZA = {theta_s}°, VZA = {theta_v}°")

        # Test with different AOD values
        for aod in [0.05, 0.1, 0.2, 0.5, 0.886]:
            r_surf = smac.smac_inv(r_toa, theta_s, phi_s, theta_v, phi_v,
                                   pressure, aod, uo3, uh2o, coefs)
            print(f"AOD={aod:.3f}: BOA reflectance = {r_surf:.4f}")

        print("\nExpected: BOA reflectance should be ~0.30-0.45 for vegetation at 850nm")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)