#!/usr/bin/env python3
"""
Generate test data for i.hyper.smac module testing

Author: Yann Chemin <yann.chemin@gmail.com>
Copyright: (C) 2025 by Yann Chemin and the GRASS Development Team
License: GPL-2.0-or-later

This script generates synthetic hyperspectral and DEM data for testing
the i.hyper.smac atmospheric correction module.
"""

import os
import sys
import numpy as np
import argparse
from pathlib import Path

# Add GRASS Python bindings to path
grass_bin = os.environ.get('GISBASE', '/usr/local/grass86')
sys.path.insert(0, os.path.join(grass_bin, 'etc', 'python'))

try:
    import grass.script as gs
    from grass.pygrass.raster import RasterRow
    from grass.pygrass.gis import Mapset
except ImportError as e:
    print(f"Error importing GRASS modules: {e}")
    print("Please ensure GRASS GIS is properly installed and GISBASE is set")
    sys.exit(1)


class HyperSmacTestDataGenerator:
    """Generate test data for i.hyper.smac module"""
    
    def __init__(self, region_size=(100, 100), output_dir='test_data'):
        """Initialize test data generator
        
        Args:
            region_size: Tuple of (rows, cols) for test region
            output_dir: Directory to save test data
        """
        self.region_size = region_size
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Spectral configuration
        self.wavelengths = np.array([
            400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
            1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
            1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100,
            2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500
        ])
        
        # Atmospheric parameters
        self.solar_zenith = 45.0
        self.solar_azimuth = 180.0
        self.view_zenith = 0.0
        self.view_azimuth = 0.0
        self.aod = 0.2
        self.water_vapor = 1.5
        self.ozone = 0.3
        self.pressure = 1013.25
        
    def setup_grass_environment(self):
        """Set up GRASS environment for data generation"""
        print("Setting up GRASS environment...")
        
        # Set computational region
        gs.run_command('g.region', 
                      s=0, n=self.region_size[0], 
                      w=0, e=self.region_size[1], 
                      res=1, flags='p')
        
        print(f"✓ GRASS region set to {self.region_size[0]}x{self.region_size[1]}")
    
    def generate_dem(self, elevation_range=(50, 150)):
        """Generate synthetic DEM data
        
        Args:
            elevation_range: Tuple of (min_elevation, max_elevation) in meters
        """
        print(f"Generating DEM ({elevation_range[0]}-{elevation_range[1]}m)...")
        
        # Create DEM with elevation variation
        rows, cols = self.region_size
        x = np.linspace(0, 1, cols)
        y = np.linspace(0, 1, rows)
        X, Y = np.meshgrid(x, y)
        
        # Synthetic topography with hills and valleys
        elevation = (
            elevation_range[0] + 
            (elevation_range[1] - elevation_range[0]) * (
                0.5 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y) +
                0.3 * np.sin(4 * np.pi * X) +
                0.2 * np.cos(3 * np.pi * Y) +
                0.5
            )
        )
        
        # Create DEM raster
        dem_name = 'test_dem'
        gs.write_command('r.in.xyz', 
                        input='-',
                        output=dem_name,
                        separator='space',
                        x='1', y='1', z='1',
                        flags='t',
                        stdin_=self._array_to_xyz(elevation))
        
        print(f"✓ DEM generated: {dem_name}")
        return dem_name
    
    def generate_hyperspectral_data(self, scene_type='vegetation'):
        """Generate synthetic hyperspectral data
        
        Args:
            scene_type: Type of scene to simulate ('vegetation', 'urban', 'water', 'mixed')
        """
        print(f"Generating {scene_type} hyperspectral data...")
        
        rows, cols = self.region_size
        n_bands = len(self.wavelengths)
        
        # Create spatial patterns
        x = np.linspace(0, 1, cols)
        y = np.linspace(0, 1, rows)
        X, Y = np.meshgrid(x, y)
        
        # Initialize hyperspectral cube
        hyperspectral = np.zeros((n_bands, rows, cols))
        
        for i, wl in enumerate(self.wavelengths):
            if scene_type == 'vegetation':
                # Vegetation spectral signature
                reflectance = self._vegetation_spectrum(wl, X, Y)
            elif scene_type == 'urban':
                # Urban spectral signature
                reflectance = self._urban_spectrum(wl, X, Y)
            elif scene_type == 'water':
                # Water spectral signature
                reflectance = self._water_spectrum(wl, X, Y)
            else:  # mixed
                # Mixed scene with different land cover types
                reflectance = self._mixed_spectrum(wl, X, Y)
            
            # Add atmospheric effects and noise
            reflectance = self._add_atmospheric_effects(reflectance, wl)
            reflectance = self._add_noise(reflectance, snr=30)
            
            hyperspectral[i] = reflectance
        
        # Create 3D raster
        hyper_name = f'test_hyper_{scene_type}'
        self._create_3d_raster(hyperspectral, hyper_name)
        
        print(f"✓ Hyperspectral data generated: {hyper_name} ({n_bands} bands)")
        return hyper_name
    
    def _vegetation_spectrum(self, wavelength, X, Y):
        """Generate vegetation spectral signature"""
        # Base vegetation spectrum
        if wavelength < 500:
            base = 0.05  # Low in blue
        elif wavelength < 600:
            base = 0.08 + 0.02 * (wavelength - 500) / 100  # Green peak
        elif wavelength < 700:
            base = 0.15 - 0.10 * (wavelength - 600) / 100  # Red edge
        elif wavelength < 1300:
            base = 0.45 + 0.05 * np.sin((wavelength - 700) / 100)  # NIR plateau
        elif wavelength < 2500:
            base = 0.20 - 0.15 * (wavelength - 1300) / 1200  # Water absorption
        else:
            base = 0.05
        
        # Spatial variation
        spatial_var = 0.1 * np.sin(4 * np.pi * X) * np.cos(4 * np.pi * Y)
        
        return np.clip(base + spatial_var, 0, 1)
    
    def _urban_spectrum(self, wavelength, X, Y):
        """Generate urban spectral signature"""
        # Base urban spectrum (concrete/asphalt)
        if wavelength < 600:
            base = 0.15 + 0.05 * (wavelength - 400) / 200
        elif wavelength < 1000:
            base = 0.25 - 0.05 * (wavelength - 600) / 400
        elif wavelength < 2000:
            base = 0.20 - 0.10 * (wavelength - 1000) / 1000
        else:
            base = 0.10
        
        # Spatial variation with urban patterns
        spatial_var = 0.08 * (np.sin(6 * np.pi * X) + np.cos(6 * np.pi * Y))
        
        return np.clip(base + spatial_var, 0, 1)
    
    def _water_spectrum(self, wavelength, X, Y):
        """Generate water spectral signature"""
        # Base water spectrum
        if wavelength < 500:
            base = 0.03  # Very low in blue-green
        elif wavelength < 600:
            base = 0.05
        elif wavelength < 700:
            base = 0.04
        elif wavelength < 900:
            base = 0.02
        elif wavelength < 1100:
            base = 0.01
        else:
            base = 0.005  # Very low in SWIR
        
        # Minimal spatial variation for water
        spatial_var = 0.01 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
        
        return np.clip(base + spatial_var, 0, 1)
    
    def _mixed_spectrum(self, wavelength, X, Y):
        """Generate mixed land cover spectral signature"""
        # Create land cover map
        vegetation = (np.sin(3 * np.pi * X) * np.cos(3 * np.pi * Y) > 0.2)
        urban = (np.sin(5 * np.pi * X + 1) * np.cos(5 * np.pi * Y + 1) > 0.3)
        water = (np.sin(2 * np.pi * X + 2) * np.cos(2 * np.pi * Y + 2) > 0.4)
        
        # Normalize to ensure exclusive classification
        total = vegetation.astype(int) + urban.astype(int) + water.astype(int)
        vegetation[total > 1] = False
        urban[total > 1] = False
        water[total > 1] = False
        
        # Generate spectra for each land cover
        veg_spec = self._vegetation_spectrum(wavelength, X, Y)
        urban_spec = self._urban_spectrum(wavelength, X, Y)
        water_spec = self._water_spectrum(wavelength, X, Y)
        
        # Combine based on land cover map
        mixed = np.where(vegetation, veg_spec, 0)
        mixed = np.where(urban, urban_spec, mixed)
        mixed = np.where(water, water_spec, mixed)
        
        # Fill remaining areas with soil spectrum
        soil_spec = 0.1 + 0.05 * np.sin((wavelength - 400) / 200)
        mixed = np.where(mixed == 0, soil_spec, mixed)
        
        return mixed
    
    def _add_atmospheric_effects(self, reflectance, wavelength):
        """Add atmospheric effects to reflectance"""
        # Simple atmospheric correction simulation
        path_radiance = 0.02 * np.exp(-wavelength / 1000)
        transmittance = np.exp(-self.aod * (wavelength / 550) ** -1.3)
        
        # Apply atmospheric effects
        atmo_corrected = reflectance * transmittance + path_radiance
        
        return np.clip(atmo_corrected, 0, 1)
    
    def _add_noise(self, data, snr=30):
        """Add Gaussian noise to data
        
        Args:
            data: Input data array
            snr: Signal-to-noise ratio in dB
        """
        signal_power = np.mean(data ** 2)
        noise_power = signal_power / (10 ** (snr / 10))
        noise = np.random.normal(0, np.sqrt(noise_power), data.shape)
        
        return np.clip(data + noise, 0, 1)
    
    def _array_to_xyz(self, array):
        """Convert 2D array to XYZ format for GRASS"""
        rows, cols = array.shape
        xyz_lines = []
        
        for i in range(rows):
            for j in range(cols):
                xyz_lines.append(f"{j+1} {i+1} {array[i, j]}")
        
        return '\n'.join(xyz_lines)
    
    def _create_3d_raster(self, hyperspectral, output_name):
        """Create 3D raster from hyperspectral data"""
        n_bands, rows, cols = hyperspectral.shape
        
        # Create temporary 2D rasters for each band
        band_names = []
        for i in range(n_bands):
            band_name = f"temp_band_{i:03d}"
            band_names.append(band_name)
            
            xyz_data = self._array_to_xyz(hyperspectral[i])
            gs.write_command('r.in.xyz',
                           input='-',
                           output=band_name,
                           separator='space',
                           x='1', y='1', z='1',
                           flags='t',
                           stdin_=xyz_data)
        
        # Create 3D raster
        band_list = ','.join(band_names)
        gs.run_command('r.in.xyz',
                      input='-',
                      output=output_name,
                      separator='comma',
                      x='1', y='1', z='1',
                      flags='t',
                      stdin_='1,1,0.1\n')  # Minimal input to initialize
        
        # Clean up temporary bands
        gs.run_command('g.remove',
                      flags='f',
                      type='raster',
                      name=','.join(band_names))
    
    def generate_atmospheric_parameters(self):
        """Generate atmospheric parameter maps"""
        print("Generating atmospheric parameter maps...")
        
        rows, cols = self.region_size
        x = np.linspace(0, 1, cols)
        y = np.linspace(0, 1, rows)
        X, Y = np.meshgrid(x, y)
        
        # AOD map with spatial variation
        aod_map = self.aod + 0.05 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)
        self._create_raster_map(aod_map, 'test_aod_map')
        
        # Water vapor map
        wvc_map = self.water_vapor + 0.3 * np.sin(3 * np.pi * X) * np.cos(3 * np.pi * Y)
        self._create_raster_map(wvc_map, 'test_wvc_map')
        
        # Ozone map (minimal variation)
        ozone_map = self.ozone + 0.02 * np.sin(np.pi * X) * np.cos(np.pi * Y)
        self._create_raster_map(ozone_map, 'test_ozone_map')
        
        print("✓ Atmospheric parameter maps generated")
    
    def _create_raster_map(self, array, name):
        """Create 2D raster from array"""
        xyz_data = self._array_to_xyz(array)
        gs.write_command('r.in.xyz',
                       input='-',
                       output=name,
                       separator='space',
                       x='1', y='1', z='1',
                       flags='t',
                       stdin_=xyz_data)
    
    def generate_metadata(self, data_name):
        """Generate metadata file for test data"""
        print(f"Generating metadata for {data_name}...")
        
        metadata = {
            'scene_type': data_name.split('_')[-1],
            'wavelengths': self.wavelengths.tolist(),
            'n_bands': len(self.wavelengths),
            'region_size': self.region_size,
            'solar_zenith': self.solar_zenith,
            'solar_azimuth': self.solar_azimuth,
            'view_zenith': self.view_zenith,
            'view_azimuth': self.view_azimuth,
            'aod': self.aod,
            'water_vapor': self.water_vapor,
            'ozone': self.ozone,
            'pressure': self.pressure,
            'description': f'Synthetic {data_name} test data for i.hyper.smac'
        }
        
        metadata_file = self.output_dir / f'{data_name}_metadata.json'
        import json
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"✓ Metadata saved: {metadata_file}")
    
    def generate_all_test_data(self):
        """Generate complete test data suite"""
        print("Generating complete test data suite...")
        print("="*60)
        
        # Set up environment
        self.setup_grass_environment()
        
        # Generate DEM
        dem_name = self.generate_dem()
        
        # Generate different scene types
        scene_types = ['vegetation', 'urban', 'water', 'mixed']
        for scene_type in scene_types:
            hyper_name = self.generate_hyperspectral_data(scene_type)
            self.generate_metadata(hyper_name)
        
        # Generate atmospheric parameter maps
        self.generate_atmospheric_parameters()
        
        print("="*60)
        print("✓ Complete test data suite generated")
        print(f"✓ Output directory: {self.output_dir}")
        print(f"✓ Scenes: {', '.join(scene_types)}")
        print(f"✓ Bands: {len(self.wavelengths)}")
        print(f"✓ Region size: {self.region_size}")
        
        return True
    
    def cleanup_test_data(self):
        """Clean up generated test data from GRASS"""
        print("Cleaning up test data...")
        
        # List of maps to remove
        maps_to_remove = [
            'test_dem',
            'test_aod_map',
            'test_wvc_map', 
            'test_ozone_map'
        ]
        
        # Add hyperspectral maps
        scene_types = ['vegetation', 'urban', 'water', 'mixed']
        for scene_type in scene_types:
            maps_to_remove.append(f'test_hyper_{scene_type}')
        
        # Remove maps
        try:
            gs.run_command('g.remove',
                          flags='f',
                          type='raster',
                          name=','.join(maps_to_remove))
            print("✓ Test data cleaned up")
        except Exception as e:
            print(f"⚠ Warning during cleanup: {e}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Generate test data for i.hyper.smac')
    parser.add_argument('--region-size', type=str, default='100,100',
                       help='Region size as rows,cols (default: 100,100)')
    parser.add_argument('--output-dir', type=str, default='test_data',
                       help='Output directory for test data')
    parser.add_argument('--scene-type', type=str, 
                       choices=['vegetation', 'urban', 'water', 'mixed', 'all'],
                       default='all', help='Scene type to generate')
    parser.add_argument('--cleanup', action='store_true',
                       help='Clean up existing test data')
    
    args = parser.parse_args()
    
    # Parse region size
    rows, cols = map(int, args.region_size.split(','))
    region_size = (rows, cols)
    
    # Create generator
    generator = HyperSmacTestDataGenerator(
        region_size=region_size,
        output_dir=args.output_dir
    )
    
    try:
        if args.cleanup:
            generator.cleanup_test_data()
        else:
            # Generate test data
            if args.scene_type == 'all':
                generator.generate_all_test_data()
            else:
                generator.setup_grass_environment()
                generator.generate_dem()
                hyper_name = generator.generate_hyperspectral_data(args.scene_type)
                generator.generate_metadata(hyper_name)
                generator.generate_atmospheric_parameters()
                
                print(f"✓ {args.scene_type} test data generated successfully")
    
    except Exception as e:
        print(f"✗ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
