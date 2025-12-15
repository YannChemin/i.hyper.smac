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
        """Find the libRadtran installation path."""
        # Check common installation paths
        paths = [
            '/usr/local/lib/libRadtran',
            '/usr/local/libRadtran',
            '/opt/libRadtran',
            '/usr/local/share/libRadtran',
            '/usr/share/libRadtran',
            os.environ.get('LIBRADTRAN_DIR', '')
        ]
        
        for path in paths:
            uvspec_path = os.path.join(path, 'bin', 'uvspec')
            if os.path.isfile(uvspec_path) and os.access(uvspec_path, os.X_OK):
                if self.verbose:
                    if self.verbose:
                        gs.message(f"Found libRadtran at: {path}")
                return path
        
        raise RuntimeError("libRadtran not found. Please install libRadtran or set LIBRADTRAN_DIR environment variable.")
    
    def _create_input_file(self, params, input_file):
        """
        Create a libRadtran input file.
        
        Args:
            params (dict): Dictionary containing the parameters for the simulation
            input_file (str): Path to the input file to create
        """
        with open(input_file, 'w') as f:
            # Atmospheric parameters
            f.write(f"atmosphere_file {os.path.join(self.libradtran_path, 'data/atmmod/afglus.dat')}\n")
            f.write(f"albedo {params.get('surface_albedo', 0.1)}\n")
            f.write(f"sza {params.get('solar_zenith', 30.0)}\n")
            f.write(f"phi0 {params.get('solar_azimuth', 0.0)}\n")
            
            # Aerosols
            if 'aerosol_model' in params:
                f.write("aerosol_default\n")
                f.write(f"aerosol_species_file {params['aerosol_model']}\n")
                f.write(f"aerosol_modify tau set {params.get('aod_550', 0.1)}\n")
            
            # Molecular absorption
            f.write("mol_abs_param reptran\n")
            f.write("mol_modify O3 {:.3f} DU\n".format(params.get('ozone', 0.3) * 1000))  # Convert to DU
            f.write("mol_modify H2O {:.3f} MOLEC\n".format(params.get('water_vapor', 2.0) * 1e-3))  # Convert to g/cm²
            
            # Wavelength grid
            f.write(f"wavelength {params['wavelength'] - params['fwhm']/2} {params['wavelength'] + params['fwhm']/2}\n")
            f.write("wavelength_grid_centers_only\n")
            
            # Output
            f.write("output_quantity transmittance\n")
            f.write("output_user lambda edir edn eup enet esum tsca tsca_dir tsca_dif\n")
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
            
        # Extract parameters
        return {
            'wavelength': float(data[0]),  # nm
            'direct_transmittance': float(data[1]),  # T_dir
            'diffuse_transmittance': float(data[2]),  # T_dif
            'spherical_albedo': float(data[3]),  # s
            'total_transmittance': float(data[4]),  # T_tot
            'path_radiance': float(data[5]),  # L_p
            'direct_irradiance': float(data[6]),  # E_dir
            'diffuse_irradiance': float(data[7]),  # E_dif
            'total_irradiance': float(data[8])  # E_tot
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
            'aerosol_model': os.path.join(
                runner.libradtran_path, 
                'data/optmod', 
                f"{aerosol_models.get(aerosol_model, 'continental_clean')}.dat"
            )
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
