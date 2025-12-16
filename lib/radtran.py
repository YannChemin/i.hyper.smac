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
from pathlib import Path

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
                        gs.message(f"Found libRadtran at: {base} (uvspec in bin/, lib/data found)")
                    # Return base path so other parts can find subdirs easily
                    return base
    
        raise RuntimeError("libRadtran not found. Ensure uvspec is in PATH and data directories exist.")

    
    def _create_input_file(self, params, input_file, transmittance=False):
        """
        Create a libRadtran input file.
        
        Args:
            params (dict): Dictionary containing the parameters for the simulation
            input_file (str): Path to the input file to create
            transmisttance (bool): False (default), if True write transmittances
        """
        with open(input_file, 'w') as f:
            dp=os.path.join(self.libradtran_path, 'share/libRadtran/data/')
            dp_atmosph=os.path.join(self.libradtran_path, 'share/libRadtran/data/atmmod/')
            #dp_solarfl=os.path.join(self.libradtran_path, 'share/libRadtran/data/solar_flux/')
            #dp_aerosol='aerosol/OPAC/standard_aerosol_files/'
            
            f.write("rte_solver sdisort\n")
            
            # Atmospheric parameters
            f.write(f"data_files_path {dp}\n")
            f.write(f"atmosphere_file {os.path.join(dp_atmosph, 'afglus.dat')}\n")
            #f.write("source solar /usr/local/share/libRadtran/data/solar_flux/kurudz_1.0nm.dat\n")
            #f.write(f"source solar {os.path.join(dp_solarfl, 'atlas_plus_modtran')}\n")
            #f.write("spline_solar_source 1\n")
            
            f.write(f"albedo {params.get('surface_albedo', 0.1)}\n")
            f.write(f"sza {params.get('solar_zenith', 30.0)}\n")
            f.write(f"phi0 {params.get('solar_azimuth', 0.0)}\n")
            
            # Aerosols
            if 'aerosol_model' in params:
                f.write("aerosol_default\n")
                f.write(f"aerosol_modify tau set {params.get('aod_550', 0.1)}\n")
            
            # Molecular absorption
            f.write("mol_abs_param reptran\n")
            f.write("mol_modify O3 {:.3f} DU\n".format(params.get('ozone', 0.3) * 1000))  # Convert to DU
            f.write("mol_modify H2O {:.3f} MM\n".format(params.get('water_vapor', 2.0)))  # Takes WVC in mm
            
            # Wavelength grid
            #low_wl= params['wavelength'] - params['fwhm']
            #hig_wl= params['wavelength'] + params['fwhm']
            #f.write(f"wavelength {low_wl} {hig_wl}\n")
            f.write(f"wavelength {params['wavelength']}\n")
            #f.write(f"wavelength_step {params['fwhm']}\n")

            # For path radiance, add: TODO satellite altitude and zenith angle
            f.write("zout 100.0\n")  # Satellite altitude in km (e.g., 100 km)
            f.write("umu 1.0\n")     # Satellite Zenith Angle (1.0 Looking straight down)
            f.write(f"phi {params.get('view_azimuth', 0.0)}\n") # Satellite Azimuth angle   

            # Output
            if transmittance == True:
                # Write default transmittance output (to be parsed)
                f.write("output_quantity transmittance\n")
            else:
                # Write Wavelength spherical albedo and path radiance 
                f.write("output_user lambda spher_alb uu\n")
            f.write("quiet\n")
    
    def run_simulation(self, params, transmittance):
        """
        Run libRadtran simulation for a single band.
        
        Args:
            params (dict): Dictionary containing the parameters for the simulation
            
        Returns:
            dict: Dictionary containing the atmospheric correction parameters
        """
        # Create input file
        input_file = os.path.join('/home/yann/RSDATA/Tanager/Kanpur/','lradtran.inp')
        self._create_input_file(params, input_file)
        input_file = os.path.join(self.temp_dir, 'lradtran.inp')
        self._create_input_file(params, input_file)
        
        # Set up command
        uvspec_path = os.path.join(self.libradtran_path, 'bin', 'uvspec')
        cmd = [uvspec_path, '<', input_file]
        
        if self.verbose:
            if self.verbose:
                gs.message(f"Running: {' '.join(cmd)}")
        
        try:
            # Run libRadtran
            result = subprocess.run(
                ' '.join(cmd),
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            
            # Parse output
            return self._parse_output(result.stdout, transmittance)
            
        except subprocess.CalledProcessError as e:
            gs.error(f"Error running libRadtran: {e}")
            gs.error(f"Stderr: {e.stderr}")
            raise
    
    def _parse_output(self, output, transmittance):
        """
        Parse libRadtran output to extract atmospheric parameters.
        
        Args:
            output (str): libRadtran output as string
            
        Returns:
            dict: Dictionary containing the atmospheric correction parameters
        """
        lines = output.strip().split('\n')
        if not lines or len(lines) < 1:
            raise ValueError("Invalid libRadtran output format")
            
        # Skip header
        data = lines[1].split()
            
        if transmittance is False:
            # Extract parameters lambda spher_alb uu
            return {
                'wavelength': float(data[0]),           # lambda (nm)
                'spherical_albedo': float(data[1]),     # spher_alb
                'path_radiance': float(data[2]),        # uu
            }
        else:
            return{
                'direct_transmittance': float(data[0]), # T_dir
                'diffuse_transmittance': float(data[1]),# T_dif
                'total_transmittance': float(data[2]),  # T_tot
        }
    
    def cleanup(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def __del__(self):
        """Destructor to ensure cleanup."""
        self.cleanup()


def get_smac_parameters(wavelength, fwhm, sza=30.0, aod_550=0.1, water_vapor=2.0, ozone=0.3, 
                       surface_albedo=0.1, aerosol_model='continental', verbose=False):
    """
    Get SMAC atmospheric correction parameters for a specific wavelength.
    
    Args:
        wavelength (float): Central wavelength in nm
        fwhm (float): Full width at half maximum in nm
        sza (float): Solar zenith angle in degrees
        aod_550 (float): Aerosol optical depth at 550nm
        water_vapor (float): Total column water vapor in g/cm²
        ozone (float): Total column ozone in cm-atm
        surface_albedo (float): Surface albedo (0-1)
        aerosol_model (str): Aerosol model (continental, maritime, urban, desert)
        verbose (bool): Enable verbose output
        
    Returns:
        dict: Dictionary containing SMAC atmospheric correction parameters
    """
    # Map aerosol model to libRadtran aerosol files
    aerosol_models = {
        'continental': 'continental_clean',
        'maritime': 'maritime_clean',
        'urban': 'urban',
        'desert': 'desert'
    }
    
    # Initialize libRadtran runner
    runner = LibRadtranRunner(verbose=verbose)
    
    try:
        # Set up parameters
        params = {
            'wavelength': wavelength,
            'fwhm': fwhm,
            'solar_zenith': sza,
            'aod_550': aod_550,
            'water_vapor': water_vapor,
            'ozone': ozone,
            'surface_albedo': surface_albedo,
            'aerosol_model': f"{aerosol_models.get(aerosol_model, 'continental_clean')}.dat"
        }
        
        # Run simulation
        transmittance = False
        result = runner.run_simulation(params, transmittance)
        transmittance = True
        resultt = runner.run_simulation(params, transmittance)
        
        # Calculate SMAC parameters
        smac_params = {
            'wavelength': wavelength,
            'direct_transmittance': resultt['direct_transmittance'],
            'diffuse_transmittance': resultt['diffuse_transmittance'],
            'spherical_albedo': result['spherical_albedo'],
            'path_radiance': result['path_radiance'],
            'total_gaseous_transmittance': resultt['direct_transmittance']  # Simplified
        }
        
        return smac_params
        
    finally:
        runner.cleanup()

def generate_smac_coefficients_libradtran(output_file, wavelengths, fwhm=10.0,
                                        solar_zenith=30.0, aod_550=0.1,
                                        water_vapor=2.0, ozone=0.3,
                                        aerosol_model='continental', verbose=True):
    """
    Generate SMAC coefficients file using libRadtran.
    
    Args:
        output_file (str): Path to save the SMAC coefficients file
        wavelengths (list): List of wavelengths (nm) to calculate coefficients for
        fwhm (float): Full width at half maximum of the spectral response function (nm)
        solar_zenith (float): Solar zenith angle in degrees
        aod_550 (float): Aerosol Optical Depth at 550nm
        water_vapor (float): Water vapor content (g/cm²)
        ozone (float): Ozone content (cm-atm)
        aerosol_model (str): Aerosol model type (e.g., 'continental', 'maritime')
        verbose (bool): Whether to print progress messages
    """
    if verbose:
        print(f"Generating SMAC coefficients using libRadtran for {len(wavelengths)} wavelengths...")
    
    # Initialize coefficient lists
    coefficients = {
        'h2o': {'a': [], 'n': []},
        'o3': {'a': [], 'n': []},
        'o2': {'a': [], 'n': [], 'p': []},
        'co2': {'a': [], 'n': [], 'p': []},
        'ch4': {'a': [], 'n': [], 'p': []},
        'no2': {'a': [], 'n': [], 'p': []},
        'co': {'a': [], 'n': [], 'p': []},
        'rayleigh': {'a0s': [], 'a1s': [], 'a2s': [], 'a3s': []},
        'aerosol': {'a0T': [], 'a1T': [], 'a2T': [], 'a3T': [], 
                   'taur': [], 'sr': [], 'a0taup': [], 'a1taup': [],
                   'wo': [], 'gc': [], 'a0P': [], 'a1P': [], 'a2P': [],
                   'a3P': [], 'a4P': [], 'Rest1': [], 'Rest2': [], 'Rest3': [],
                   'Rest4': [], 'Resr1': [], 'Resr2': [], 'Resr3': [],
                   'Resa1': [], 'Resa2': [], 'Resa3': [], 'Resa4': []}
    }
    
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        
        # Run libRadtran for each wavelength
        for wl in wavelengths:
            if verbose:
                print(f"Processing wavelength: {wl} nm")
            
            # Create input file for libRadtran
            input_file = tmp_dir / f"libradtran_input_{wl}.inp"
            output_file_tmp = tmp_dir / f"libradtran_output_{wl}.out"
            
            # Write libRadtran input file
            with open(input_file, 'w') as f:
                f.write(f"wavelength {wl-5} {wl+5}\n")  # Small range around the wavelength
                f.write(f"solar_file data/solar_flux/kurudz_0.1nm.dat\n")
                f.write(f"atmosphere_file data/atmmod/afglms.dat\n")
                f.write(f"source solar\n")
                f.write(f"day_of_year 180\n")
                f.write(f"latitude 0.0\n")
                f.write(f"longitude 0.0\n")
                f.write(f"time 12:00:00\n")
                f.write(f"sza {solar_zenith}\n")
                f.write(f"phi 0\n")
                f.write(f"umu 1.0\n")
                f.write(f"phi0 0\n")
                f.write(f"altitude 0.0\n")
                f.write(f"mol_abs_param reptran medium\n")
                f.write(f"rte_solver disort\n")
                f.write(f"number_of_streams 8\n")
                f.write(f"aerosol_default\n")
                f.write(f"aerosol_angstrom 1.3 0.0\n")
                f.write(f"aerosol_modify tau set {aod_550} 550\n")
                f.write(f"aerosol_haze 1\n")
                f.write(f"aerosol_season 1\n")
                f.write(f"aerosol_angstrom_exponent 1.3\n")
                f.write(f"aerosol_angstrom_lambda1 400.0\n")
                f.write(f"aerosol_angstrom_lambda2 900.0\n")
                f.write(f"mol_modify H2O {water_vapor} MM\n")
                f.write(f"mol_modify O3 {ozone} DU\n")
                f.write(f"output_user lambda eglo\n")
                f.write(f"output_quantity transmittance\n")
                f.write(f"quiet\n")
            
            # Run libRadtran
            try:
                cmd = f"uvspec < {input_file} > {output_file_tmp}"
                subprocess.run(cmd, shell=True, check=True)
                
                # Parse output and extract coefficients
                # This is a simplified example - you'll need to parse the actual output
                # and map it to SMAC coefficients
                with open(output_file_tmp, 'r') as f:
                    lines = f.readlines()
                    # Parse the output and extract necessary values
                    # This is a placeholder - actual parsing needed
                    pass
                    
                # Example of how you might extract coefficients
                # (replace with actual parsing logic)
                coefficients['h2o']['a'].append(0.0)  # Replace with actual value
                coefficients['h2o']['n'].append(0.0)
                # ... (add other coefficients)
                
            except subprocess.CalledProcessError as e:
                print(f"Error running libRadtran for wavelength {wl} nm: {e}")
                continue
    
    # After processing all wavelengths, calculate mean values and write to file
    with open(output_file, 'w') as f:
        # H2O line
        f.write(f"{np.mean(coefficients['h2o']['a']):.6f} {np.mean(coefficients['h2o']['n']):.6f}\n")
        
        # O3 line
        f.write(f"{np.mean(coefficients['o3']['a']):.6f} {np.mean(coefficients['o3']['n']):.6f}\n")
        
        # O2 line
        f.write(f"{np.mean(coefficients['o2']['a']):.6f} {np.mean(coefficients['o2']['n']):.6f} {np.mean(coefficients['o2']['p']):.6f}\n")
        
        # CO2 line
        f.write(f"{np.mean(coefficients['co2']['a']):.6f} {np.mean(coefficients['co2']['n']):.6f} {np.mean(coefficients['co2']['p']):.6f}\n")
        
        # CH4 line
        f.write(f"{np.mean(coefficients['ch4']['a']):.6f} {np.mean(coefficients['ch4']['n']):.6f} {np.mean(coefficients['ch4']['p']):.6f}\n")
        
        # NO2 line
        f.write(f"{np.mean(coefficients['no2']['a']):.6f} {np.mean(coefficients['no2']['n']):.6f} {np.mean(coefficients['no2']['p']):.6f}\n")
        
        # CO line
        f.write(f"{np.mean(coefficients['co']['a']):.6f} {np.mean(coefficients['co']['n']):.6f} {np.mean(coefficients['co']['p']):.6f}\n")
        
        # Rayleigh and aerosol scattering (lines 8-19)
        # These are placeholders - replace with actual calculations
        f.write(f"{0.0:.6f} {0.0:.6f} {0.0:.6f} {0.0:.6f}\n")  # Line 8
        f.write(f"{0.0:.6f} {0.0:.6f} {0.0:.6f} {0.0:.6f}\n")  # Line 9
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 10
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 11
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 12
        f.write(f"{0.0:.6f} {0.0:.6f} {0.0:.6f}\n")  # Line 13
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 14
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 15
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 16
        f.write(f"{0.0:.6f} {0.0:.6f} {0.0:.6f}\n")  # Line 17
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 18
        f.write(f"{0.0:.6f} {0.0:.6f}\n")  # Line 19
    
    if verbose:
        print(f"SMAC coefficients written to: {output_file}")
    
    return output_file

if __name__ == "__main__":
    # Example usage
    params = get_smac_parameters(
        wavelength=550.0,
        fwhm=10.0,
        sza=30.0,
        aod_550=0.2,
        water_vapor=2.5,
        ozone=0.3,
        surface_albedo=0.1,
        aerosol_model='continental',
        verbose=True
    )
    if gs.verbosity() > 0:
        gs.message(f"SMAC Parameters: {params}")
