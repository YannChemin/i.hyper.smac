# i.hyper.smac Programmer Manual

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Architecture](#architecture)
4. [API Reference](#api-reference)
5. [Testing](#testing)
6. [Performance Optimization](#performance-optimization)
7. [Troubleshooting](#troubleshooting)
8. [Contributing](#contributing)

## Overview

The `i.hyper.smac` module provides advanced atmospheric correction capabilities for hyperspectral imagery using the Simplified Method for Atmospheric Correction (SMAC) enhanced with libRadtran DISORT calculations. This module implements state-of-the-art atmospheric correction algorithms including automatic aerosol and water vapor estimation, parallel processing, and OpenCL acceleration.

### Key Features

- **Full libRadtran DISORT Integration**: Multiple scattering atmospheric correction
- **Automatic AOD Estimation**: Dark Target retrieval algorithm
- **Water Vapor Retrieval**: Joint optimal estimation and band ratio methods
- **Smart LUT System**: Scene-specific Look-Up Table optimization (77% speedup)
- **Parallel Processing**: Multi-core LUT generation with automatic resource detection
- **OpenCL Acceleration**: GPU/CPU acceleration for atmospheric correction
- **Uncertainty Quantification**: Per-band uncertainty propagation
- **Adjacency Correction**: Point Spread Function based correction
- **Spectral Polishing**: ISOFIT-inspired spectral smoothing

## Installation

### System Requirements

- **GRASS GIS**: Version 8.0 or higher
- **Python**: 3.7 or higher
- **libRadtran**: For atmospheric modeling
- **OpenCL**: For GPU acceleration (optional)

### Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install grass grass-dev libradtran-dev libopencl-dev

# Python dependencies
pip install numpy scipy pyopencl netcdf4 h5py
```

### Module Installation

```bash
# Clone repository
git clone https://github.com/YannChemin/i.hyper.smac.git
cd i.hyper.smac

# Install module
make install

# Test installation
make test
```

## Architecture

### System Architecture

```
i.hyper.smac (main module)
├── Atmospheric Parameter Estimation
│   ├── Dark Target AOD Retrieval
│   ├── Water Vapor Content Retrieval
│   └── Pressure Estimation from DEM
├── LUT Generation and Management
│   ├── Smart LUT System (smart_lut.py)
│   ├── Parallel LUT Generator (parallel_lut.py)
│   └── Standard LUT (lut.py)
├── Atmospheric Correction Engine
│   ├── OpenCL Accelerator (opencl_accelerator.py)
│   ├── Spectral Polishing
│   └── Adjacency Correction
└── Uncertainty Quantification
    ├── Error Propagation
    └── MAP Inversion
```

### Data Flow

1. **Input Processing**: Read hyperspectral 3D raster and metadata
2. **Atmospheric Estimation**: Estimate AOD, water vapor, and pressure
3. **LUT Generation**: Create or load optimized Look-Up Tables
4. **Atmospheric Correction**: Apply correction using OpenCL acceleration
5. **Post-Processing**: Apply spectral polishing and adjacency correction
6. **Output Generation**: Write corrected data and uncertainty maps

## API Reference

### Main Functions

#### `apply_lut_correction()`

Main atmospheric correction function.

```python
def apply_lut_correction(input_raster, output_raster, bands,
                         aod, water_vapor, ozone, pressure,
                         solar_zenith, solar_azimuth,
                         view_zenith, view_azimuth,
                         aerosol_model='continental', 
                         smart_lut='auto',
                         parallel_lut='auto',
                         opencl_device='auto',
                         **kwargs):
    """
    Apply atmospheric correction to hyperspectral data
    
    Parameters:
    -----------
    input_raster : str
        Input hyperspectral 3D raster name
    output_raster : str  
        Output corrected 3D raster name
    bands : list
        List of band wavelength dictionaries
    aod : float
        Aerosol Optical Depth at 550nm
    water_vapor : float or str
        Water vapor content (g/cm²) or retrieval method
    ozone : float
        Ozone content (cm-atm)
    pressure : float
        Atmospheric pressure (hPa)
    solar_zenith : float
        Solar zenith angle (degrees)
    solar_azimuth : float
        Solar azimuth angle (degrees)
    view_zenith : float
        View zenith angle (degrees)
    view_azimuth : float
        View azimuth angle (degrees)
    aerosol_model : str
        Aerosol model type
    smart_lut : str
        Smart LUT mode ('auto', 'yes', 'no')
    parallel_lut : str
        Parallel LUT mode ('auto', 'enabled', 'disabled')
    opencl_device : str
        OpenCL device type ('auto', 'gpu', 'cpu')
    
    Returns:
    --------
    AtmosphericLUT object
        Generated or loaded atmospheric LUT
    """
```

### LUT Classes

#### `AtmosphericLUT`

Atmospheric Look-Up Table for hyperspectral correction.

```python
class AtmosphericLUT:
    """
    Atmospheric Look-Up Table for hyperspectral correction
    
    Attributes:
    -----------
    wavelengths : array
        Wavelength array (nm)
    aod_grid : array
        Aerosol Optical Depth grid
    h2o_grid : array
        Water vapor grid (g/cm²)
    r_atm : array
        Atmospheric reflectance (wavelength × AOD × H2O)
    t_down : array
        Downward transmittance
    t_up : array
        Upward transmittance
    s : array
        Spherical albedo
    """
    
    @classmethod
    def get_or_generate(cls, sza, vza, phi, pressure, aerosol_model,
                       wl_min, wl_max, wl_step, h2o=None, o3=None,
                       **kwargs):
        """
        Get cached LUT or generate new one
        
        Implements intelligent caching with scene-aware keys
        """
```

#### `SmartLUTManager`

Smart LUT management system.

```python
class SmartLUTManager:
    """
    Smart LUT management system
    
    Methods:
    --------
    get_smart_lut_config(scene_aod, precision_threshold)
        Generate scene-specific AOD grid configuration
    
    should_regenerate_lut(scene_aod, cache_path, precision_threshold)
        Check if LUT regeneration is needed
    
    generate_cache_key(scene_parameters)
        Generate cache key for LUT storage
    """
```

## Testing

### Test Suite Structure

```
testsuite/
├── test_i_hyper_smac.py    # Main unit test file
├── run_tests.py            # Comprehensive test runner
├── generate_test_data.py   # Test data generator
├── Makefile               # GRASS make integration
└── README.md              # Test documentation
```

### Running Tests

```bash
# Run all tests
make test-all

# Run specific test categories
make test-unit          # Unit tests only
make test-integration    # Integration tests only
make test-performance    # Performance tests only

# Generate test coverage
make test-coverage

# Run with custom test data
cd testsuite
python run_tests.py --test-data /path/to/test/data
```

### Test Categories

#### Unit Tests

- **Basic Functionality**: Module existence and basic operations
- **Parameter Validation**: Input validation and error handling
- **AOD Estimation**: Dark Target algorithm testing
- **Water Vapor Retrieval**: Multiple retrieval methods
- **Smart LUT System**: Performance optimization validation
- **OpenCL Acceleration**: GPU/CPU acceleration testing

#### Integration Tests

- **End-to-End Workflows**: Complete atmospheric correction pipelines
- **Output Generation**: File creation and verification
- **Error Handling**: Graceful failure and recovery

#### Performance Tests

- **Smart LUT vs Standard LUT**: ~77% speedup validation
- **Parallel vs Serial Processing**: Multi-core acceleration
- **OpenCL vs CPU**: GPU acceleration benchmarks

## Performance Optimization

### Smart LUT System

The Smart LUT system optimizes LUT generation based on scene characteristics:

```python
# Standard LUT: 96 runs (13 AOD × 8 H2O × 2 albedo)
# Smart LUT: 36 runs (3-5 AOD × 6 H2O × 2 albedo)
# Speedup: ~77% reduction in computation time
```

### Parallel Processing

Multi-core LUT generation with automatic resource detection:

```python
# Enable parallel processing
gs.run_command('i.hyper.smac',
               input='hyperspectral',
               output='corrected',
               parallel_lut='enabled',
               max_workers='auto')
```

### OpenCL Acceleration

GPU/CPU acceleration for atmospheric correction:

```python
# GPU acceleration
gs.run_command('i.hyper.smac',
               input='hyperspectral',
               output='corrected',
               opencl_device='gpu',
               opencl_memory=2048)
```

## Troubleshooting

### Common Issues

#### libRadtran not found

```bash
# Install libRadtran
sudo apt-get install libradtran-dev

# Set environment variable
export LIBRADTRAN_PATH=/usr/local/lib
```

#### OpenCL initialization fails

```bash
# Install OpenCL drivers
sudo apt-get install opencl-headers ocl-icd-opencl-dev

# Use CPU fallback
gs.run_command('i.hyper.smac', opencl_device='cpu')
```

#### Memory errors

```bash
# Reduce memory usage
gs.run_command('i.hyper.smac', opencl_memory=512)
```

### Debug Mode

```bash
# Enable verbose output
export GRASS_VERBOSE=3

# Run with debugging
gs.run_command('i.hyper.smac', flags='v')
```

## Contributing

### Development Setup

```bash
# Clone repository
git clone https://github.com/YannChemin/i.hyper.smac.git
cd i.hyper.smac

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
make test-all
```

### Code Style

Follow GRASS GIS coding standards:

- **Python**: PEP 8 with 4-space indentation
- **Documentation**: Doxygen format for C/C++, docstrings for Python
- **Testing**: Comprehensive test coverage required

### Submitting Changes

1. Fork repository
2. Create feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit pull request

### Test Coverage

Maintain high test coverage:

```bash
# Generate coverage report
make test-coverage

# View coverage report
open testsuite/htmlcov/index.html
```

## References

- Vermote, E., et al. (1997). Second simulation of the satellite signal in the solar spectrum (6S). IEEE TGRS.
- Levy, R.C., et al. (2007). Global aerosol optical properties and application to MODIS aerosol retrieval over land. JGR.
- Thompson, D.R., et al. (2018). ISOFIT: A full-physics, optimal estimation retrieval package for remote sensing. JQSRT.
- Mayer, B., & Kylling, A. (2005). Technical note: The libRadtran software package for radiative transfer calculations. ACP.

## License

This program is free software under GNU General Public License (>=v2). 
Read the file COPYING that comes with GRASS for details.

## Support

- **GRASS GIS Website**: https://grass.osgeo.org
- **Documentation**: https://grass.osgeo.org/documentation/
- **Mailing List**: grass-dev@lists.osgeo.org
- **Issue Tracker**: https://github.com/YannChemin/i.hyper.smac/issues
