#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
libRadtran interface for atmospheric correction parameter generation.

This module provides functionality to generate atmospheric correction parameters
using libRadtran for SMAC atmospheric correction.
"""

import os
import tempfile
import subprocess
import grass.script as gs
import numpy as np
import sys
from datetime import datetime

# Import smac module for coefficient class
try:
    import smac
except ImportError:
    gs.warning("Could not import smac module directly, will attempt relative import")
    from . import smac

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
                f.write(f"aerosol_modify tau set {params.get('aod_550', 0.1)}\n")
            
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
            return {
                'direct_transmittance': float(data[0]),
                'diffuse_transmittance': float(data[1]),
                'total_transmittance': float(data[2]),
            }
        elif output_type == 'radiance':
            return {
                'wavelength': float(data[0]),
                'spherical_albedo': float(data[1]),
                'path_radiance': float(data[2]),
            }
        else:  # reflectance
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
       atmosphere_file="/usr/local/share/libRadtran/data/atmmod/afglus.dat"):
    """
    Compute band-integrated exo-atmospheric irradiance E0_band using libRadtran.
    """
    # Calculate wavelength range (nm)
    wl_min = wl_center - 2 * fwhm
    wl_max = wl_center + 2 * fwhm

    # Create a temporary directory
    import tempfile
    import os
    import numpy as np

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, "E0_spectrum.dat")
        
        # Build uvspec input file (TOA solar irradiance spectrum)
        uvspec_inp = f"""\
atmosphere_file {atmosphere_file}
source solar {solar_file} per_nm
wavelength {wl_min} {wl_max}
spline
output_quantity irradiance
output_user lambda e_dir
output_process sum
output process_quantity per_nm
output_skip_irradiance_scaling 1
verbose
"""
        # Write input file
        inp_path = os.path.join(tmpdir, "uvspec.inp")
        with open(inp_path, 'w') as f:
            f.write(uvspec_inp)

        # Run uvspec
        import subprocess
        try:
            result = subprocess.run(
                [uvspec_bin, "<", inp_path],
                capture_output=True, text=True, shell=True
            )
            if result.returncode != 0:
                print(f"Error running uvspec: {result.stderr}")
                return None
        except Exception as e:
            print(f"Error: {e}")
            return None

        # Process output
        data = np.loadtxt(out_path)
        lam_uv = data[:, 0]    # nm
        edir = data[:, 1]      # direct solar irradiance at TOA

    # Interpolate to regular wavelength grid, if needed
    wavelengths = np.linspace(lam_uv.min(), lam_uv.max(), 1000)
    E0_lambda = np.interp(wavelengths, lam_uv, edir)

    # Build band response and integrate
    R = gaussian_rsp(wavelengths, wl_center, fwhm)
    num = np.trapz(E0_lambda * R, wavelengths)
    den = np.trapz(R, wavelengths)
    E0_band = num / den

    return E0_band

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
    Generate SMAC-compatible coefficient object from libRadtran simulations.
    
    This function runs multiple libRadtran simulations to estimate the coefficients
    needed by the SMAC atmospheric correction model.
    
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
        smac.coeff: SMAC coefficient object compatible with smac.py
    """
    
    if verbose:
        gs.message(f"Generating SMAC coefficients for {wavelength} nm using libRadtran...")
    
    runner = LibRadtranRunner(verbose=verbose)
    
    try:
        # Base parameters for all simulations
        base_params = {
            'wavelength': wavelength,
            'solar_zenith': solar_zenith,
            'solar_azimuth': solar_azimuth,
            'view_azimuth': view_azimuth,
            'aod_550': aod_550,
            'water_vapor': water_vapor,
            'ozone': ozone,
            'pressure': pressure,
            'surface_albedo': 0.1,
            'aerosol_model': aerosol_model
        }
        
        # Run simulations to get transmittances
        trans_result = runner.run_simulation(base_params, 'transmittance')
        
        # Run simulation for spherical albedo and path radiance
        rad_result = runner.run_simulation(base_params, 'radiance')
        
        # Create a temporary SMAC coefficient file
        temp_coef_file = os.path.join(runner.temp_dir, f'smac_coef_{wavelength}.dat')
        
        # Generate approximate SMAC coefficients
        # These are simplified approximations based on libRadtran output
        with open(temp_coef_file, 'w') as f:
            # H2O coefficients (line 1)
            ah2o = -0.1 if 850 < wavelength < 1050 else -0.01
            f.write(f"{ah2o:.6f} 0.5\n")
            
            # O3 coefficients (line 2)
            ao3 = -0.05 if 400 < wavelength < 700 else 0.0
            f.write(f"{ao3:.6f} 0.5\n")
            
            # O2 coefficients (line 3)
            f.write("0.0 0.0 1.0\n")
            
            # CO2 coefficients (line 4)
            f.write("0.0 0.0 1.0\n")
            
            # CH4 coefficients (line 5)
            f.write("0.0 0.0 1.0\n")
            
            # NO2 coefficients (line 6)
            f.write("0.0 0.0 1.0\n")
            
            # CO coefficients (line 7)
            f.write("0.0 0.0 1.0\n")
            
            # Rayleigh scattering coefficients (line 8)
            # s = a0s*P + a3s + a1s*tau + a2s*tau^2
            a0s = rad_result['spherical_albedo'] / (pressure / 1013.25)
            a1s = 0.0
            a2s = 0.0
            a3s = 0.0
            f.write(f"{a0s:.6f} {a1s:.6f} {a2s:.6f} {a3s:.6f}\n")
            
            # Transmission coefficients (line 9)
            # T = a0T + a1T*tau/cos(theta) + (a2T*P + a3T)/(1+cos(theta))
            a0T = trans_result['direct_transmittance']
            a1T = 0.0
            a2T = 0.0
            a3T = 0.0
            f.write(f"{a0T:.6f} {a1T:.6f} {a2T:.6f} {a3T:.6f}\n")
            
            # Rayleigh optical thickness (line 10)
            taur = 0.008569 * (wavelength / 1000) ** (-4)
            f.write(f"{taur:.6f} {taur:.6f}\n")
            
            # Aerosol optical thickness coefficients (line 11)
            # taup = a0taup + a1taup*taup550
            alpha = 1.3  # Angstrom exponent
            a0taup = 0.0
            a1taup = (wavelength / 550.0) ** (-alpha)
            f.write(f"{a0taup:.6f} {a1taup:.6f}\n")
            
            # Aerosol single scattering albedo and asymmetry parameter (line 12)
            wo = 0.9  # Typical single scattering albedo
            gc = 0.7  # Typical asymmetry parameter
            f.write(f"{wo:.6f} {gc:.6f}\n")
            
            # Phase function coefficients (line 13-14)
            f.write("0.0 0.0 0.0\n")
            f.write("0.0 0.0\n")
            
            # Residual terms (lines 15-19)
            # These account for coupling effects and residuals
            f.write("0.0 0.0\n")  # Rest1, Rest2
            f.write("0.0 0.0\n")  # Rest3, Rest4
            f.write("0.0 0.0 0.0\n")  # Resr1, Resr2, Resr3
            f.write("0.0 0.0\n")  # Resa1, Resa2
            f.write("0.0 0.0\n")  # Resa3, Resa4
        
        # Load the coefficients using smac.coeff class
        coefs = smac.coeff(temp_coef_file)
        
        return coefs
        
    finally:
        runner.cleanup()


def get_smac_parameters(wavelength, fwhm=10.0, sza=30.0, vza=0.0, 
                       aod_550=0.1, water_vapor=2.0, ozone=0.3,
                       pressure=1013.25, surface_albedo=0.1, 
                       aerosol_model='continental', verbose=False):
    """
    Get SMAC atmospheric correction parameters for a specific wavelength.
    
    This function generates a full set of SMAC coefficients and returns them
    in a format compatible with the smac_inv function.
    
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
        verbose (bool): Enable verbose output
        
    Returns:
        smac.coeff: SMAC coefficient object for use with smac_inv
    """
    
    coefs = generate_smac_coefficients_from_libradtran(
        wavelength=wavelength,
        fwhm=fwhm,
        solar_zenith=sza,
        view_zenith=vza,
        aod_550=aod_550,
        water_vapor=water_vapor,
        ozone=ozone,
        pressure=pressure,
        aerosol_model=aerosol_model,
        verbose=verbose
    )
    
    return coefs


if __name__ == "__main__":
    # Example usage
    import sys
    
    # Test the coefficient generation
    wavelength = 550.0
    
    try:
        coefs = get_smac_parameters(
            wavelength=wavelength,
            fwhm=10.0,
            sza=30.0,
            aod_550=0.2,
            water_vapor=2.5,
            ozone=0.3,
            surface_albedo=0.1,
            aerosol_model='continental',
            verbose=True
        )
        
        print(f"\nGenerated SMAC coefficients for {wavelength} nm")
        print(f"H2O: ah2o={coefs.ah2o:.6f}, nh2o={coefs.nh2o:.6f}")
        print(f"O3: ao3={coefs.ao3:.6f}, no3={coefs.no3:.6f}")
        print(f"Rayleigh: taur={coefs.taur:.6f}")
        print(f"Aerosol: wo={coefs.wo:.6f}, gc={coefs.gc:.6f}")
        
        # Test the smac_inv function
        r_toa = 0.2
        theta_s = 30.0
        phi_s = 0.0
        theta_v = 0.0
        phi_v = 0.0
        pressure = 1013.25
        taup550 = 0.1
        uo3 = 0.3
        uh2o = 2.0
        
        r_surf = smac.smac_inv(r_toa, theta_s, phi_s, theta_v, phi_v,
                               pressure, taup550, uo3, uh2o, coefs)
        
        print(f"\nTest atmospheric correction:")
        print(f"TOA reflectance: {r_toa}")
        print(f"Surface reflectance: {r_surf:.4f}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)