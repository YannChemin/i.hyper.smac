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
            # Atmospheric parameters
            f.write(f"data_files_path {os.path.join(self.libradtran_path, 'share/libRadtran/data/')}\n")
            f.write(f"atmosphere_file {os.path.join(self.libradtran_path, 'share/libRadtran/data/atmmod/afglus.dat')}\n")
            f.write(f"source solar {os.path.join(self.libradtran_path, 'share/libRadtran/data/solar_flux/atlas_plus_modtran')}\n")
            
            f.write(f"albedo {params.get('surface_albedo', 0.1)}\n")
            f.write(f"sza {params.get('solar_zenith', 30.0)}\n")
            f.write(f"phi0 {params.get('solar_azimuth', 0.0)}\n")
            
            # Aerosols
            if 'aerosol_model' in params:
                f.write("aerosol_default\n")
                f.write(f"aerosol_species_file {os.path.join(self.libradtran_path, 'share/libRadtran/data/', params['aerosol_model'])}\n")
                f.write(f"aerosol_modify tau set {params.get('aod_550', 0.1)}\n")
            
            # Molecular absorption
            f.write("mol_abs_param reptran\n")
            f.write("mol_modify O3 {:.3f} DU\n".format(params.get('ozone', 0.3) * 1000))  # Convert to DU
            f.write("mol_modify H2O {:.3f} MM\n".format(params.get('water_vapor', 2.0)))  # Takes WVC in mm
            
            # Wavelength grid
            low_wl= params['wavelength'] - params['fwhm']
            hig_wl= params['wavelength'] + params['fwhm']
            f.write(f"wavelength {low_wl} {hig_wl}\n")
            #f.write(f"wavelength_step {params['fwhm']}\n")

            # For path radiance, add: TODO satellite altitude and zenith angle
            f.write("zout 100.0\n")  # Satellite altitude in km (e.g., 100 km)
            f.write("umu 1.0\n")     # Satellite Zenith Angle (1.0 Looking straight down)
            f.write(f"phi {params.get('view_azimuth', 0.0)}\n") # Satellite Azimuth angle   

            # Output
            if transmittance == True:
                f.write("output_quantity transmittance\n")
            else:
                f.write("output_user lambda edir edn eup enet eglo spher_alb uu\n")
            f.write("quiet\n")
    
    def run_simulation(self, params):
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
            return self._parse_output(result.stdout)
            
        except subprocess.CalledProcessError as e:
            gs.error(f"Error running libRadtran: {e}")
            gs.error(f"Stderr: {e.stderr}")
            raise
    
    def _parse_output(self, output):
        """
        Parse libRadtran output to extract atmospheric parameters.
        
        Args:
            output (str): libRadtran output as string
            
        Returns:
            dict: Dictionary containing the atmospheric correction parameters
        """
        lines = output.strip().split('\n')
        if not lines or len(lines) < 2:
            raise ValueError("Invalid libRadtran output format")
            
        # Skip header
        data = lines[1].split()
        if len(data) < 9:
            raise ValueError("Unexpected number of columns in libRadtran output")
            
        # Extract parameters lambda edir edn eup enet eglo spher_alb uu
        return {
            'wavelength': float(data[0]),           # lambda (nm)
            'direct_irradiance': float(data[1]),    # edir
            'diffuse_irradiance': float(data[2]),   # edn (downward)
            'upward_irradiance': float(data[3]),    # eup
            'net_irradiance': float(data[4]),       # enet
            'global_irradiance': float(data[5]),    # eglo = edir + edn
            'spherical_albedo': float(data[6]),     # spher_alb
            'path_radiance': float(data[7]),        # uu

            'direct_transmittance': float(data[8]),  # T_dir
            'diffuse_transmittance': float(data[9]),  # T_dif
            'total_transmittance': float(data[10]),  # T_tot
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
        result = runner.run_simulation(params)
        
        # Calculate SMAC parameters
        smac_params = {
            'wavelength': wavelength,
            'direct_transmittance': result['direct_transmittance'],
            'diffuse_transmittance': result['diffuse_transmittance'],
            'spherical_albedo': result['spherical_albedo'],
            'path_radiance': result['path_radiance'],
            'total_gaseous_transmittance': result['direct_transmittance']  # Simplified
        }
        
        return smac_params
        
    finally:
        runner.cleanup()


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
