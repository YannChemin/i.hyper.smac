# i.hyper.smac Test Suite

This directory contains comprehensive tests for the i.hyper.smac hyperspectral atmospheric correction module, following GRASS GIS testing standards.

## Overview

The test suite provides:
- **Unit Tests**: Individual function and component testing
- **Integration Tests**: End-to-end workflow testing
- **Performance Tests**: Benchmarking and optimization validation
- **Programmer Manual**: Complete API documentation

## Test Structure

```
testsuite/
├── test_i_hyper_smac.py    # Main unit test file
├── run_tests.py            # Comprehensive test runner
├── Makefile               # GRASS make integration
├── README.md              # This file
└── test_data/             # Test data directory (optional)
    ├── sample_hyper.tif   # Sample hyperspectral data
    └── sample_dem.tif     # Sample DEM data
```

## Running Tests

### Quick Test Run

```bash
# Run all tests
python run_tests.py

# Run only unit tests
python run_tests.py --unit-only

# Run with verbose output
python run_tests.py --verbose
```

### GRASS Integration

```bash
# Run using GRASS make system
make parsertest

# Run individual test file
grass --tmp-location EPSG:4326 --exec python test_i_hyper_smac.py
```

### Test Categories

#### 1. Unit Tests (`test_i_hyper_smac.py`)

**TestHyperSmacBasic**
- Module existence and basic functionality
- Required parameter validation
- Basic atmospheric correction
- Parameter range validation
- Smart LUT functionality
- Parallel LUT functionality

**TestHyperSmacAODEstimation**
- Automatic AOD estimation
- AOD map generation
- Dark Target algorithm validation

**TestHyperSmacWaterVapor**
- Water vapor estimation methods
- WVC map generation
- Joint optimal estimation

**TestHyperSmacAdvanced**
- Uncertainty estimation
- Adjacency correction
- Aerosol models
- OpenCL functionality
- Cache regeneration

**TestHyperSmacPerformance**
- Smart vs Standard LUT performance
- Parallel vs Serial LUT performance
- Resource usage optimization

#### 2. Integration Tests

- **Basic Integration**: End-to-end workflow testing
- **Parameter Validation**: Input validation and error handling
- **Output Generation**: Output file creation and verification

#### 3. Performance Tests

- **Smart LUT Performance**: ~77% speedup validation
- **Parallel LUT Performance**: Multi-core acceleration
- **OpenCL Performance**: GPU/CPU acceleration

## Test Requirements

### System Requirements

- **GRASS GIS**: Version 8.0 or higher
- **Python**: 3.7 or higher
- **libRadtran**: For atmospheric modeling
- **OpenCL**: For GPU acceleration (optional)

### Python Dependencies

```bash
# Required GRASS Python bindings
grass.script
grass.gunittest

# Additional dependencies
numpy
scipy
pyopencl  # For OpenCL tests
```

### Test Data

Tests create synthetic data automatically, but real data can be used:

```bash
# Provide custom test data directory
python run_tests.py --test-data /path/to/test/data
```

## Test Configuration

### Environment Variables

```bash
export GRASS_VERBOSE=3          # Maximum verbosity
export LIBRADTRAN_PATH=/usr/local/lib
export OPENCL_DEVICE_TYPE=gpu   # Preferred OpenCL device
export SMART_LUT_CACHE_DIR=/tmp/smart_lut_cache
```

### Performance Tuning

```python
# In test_i_hyper_smac.py
PERFORMANCE_CONFIG = {
    'opencl_device': 'gpu',
    'opencl_memory': 2048,
    'parallel_lut': 'enabled',
    'smart_lut': 'yes',
    'max_workers': 'auto',
    'resource_usage': 0.75
}
```

## Expected Results

### Unit Tests

- **Total Tests**: ~25 test methods
- **Execution Time**: 2-5 minutes
- **Coverage**: Core functionality and edge cases

### Integration Tests

- **Test Scenarios**: 5-8 workflows
- **Execution Time**: 1-3 minutes
- **Data Size**: Small synthetic datasets

### Performance Tests

- **Smart LUT**: Should show ~77% speedup
- **Parallel LUT**: Should utilize multiple cores
- **OpenCL**: Should show GPU acceleration (if available)

## Troubleshooting

### Common Issues

1. **libRadtran not found**
   ```bash
   sudo apt-get install libradtran-dev
   export LIBRADTRAN_PATH=/usr/local/lib
   ```

2. **OpenCL initialization fails**
   ```bash
   # Install OpenCL drivers
   sudo apt-get install opencl-headers ocl-icd-opencl-dev
   
   # Test with CPU fallback
   python run_tests.py --verbose
   ```

3. **Memory errors**
   ```bash
   # Reduce memory usage
   export OPENCL_MEMORY=512
   ```

4. **Permission errors**
   ```bash
   # Ensure write permissions for test directory
   chmod 755 testsuite/
   ```

### Debug Mode

```bash
# Enable maximum verbosity
export GRASS_VERBOSE=3

# Run with Python debugger
python -m pdb run_tests.py --unit-only
```

### Test Data Issues

```bash
# Check test data integrity
python -c "
import grass.script as gs
gs.run_command('r.in.xyz', input='test_data.csv', output='test', flags='t')
print('Test data loaded successfully')
"
```

## Contributing Tests

### Adding New Tests

1. **Create test class** inheriting from `TestCase`
2. **Follow naming convention**: `TestFeatureName`
3. **Use descriptive test methods**: `test_specific_functionality`
4. **Include setup/teardown**: `setUpClass()`, `tearDownClass()`
5. **Use GRASS assertions**: `assertModule()`, `assertRasterExists()`

### Example Test Structure

```python
class TestNewFeature(TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment"""
        cls.runModule("g.region", s=0, n=10, w=0, e=10, res=1)
        cls.use_temp_region()
        # Create test data...

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment"""
        # Remove test maps...
        cls.del_temp_region()

    def test_new_functionality(self):
        """Test new feature"""
        self.assertModule("i.hyper.smac", new_parameter="value")
        self.assertTrue(is_map_in_mapset("expected_output"))
```

### Test Guidelines

- **Isolation**: Tests should not depend on each other
- **Cleanup**: Always clean up temporary data
- **Assertions**: Use appropriate GRASS assertions
- **Documentation**: Include docstrings for test methods
- **Performance**: Keep test execution time reasonable

## Continuous Integration

### GitHub Actions Example

```yaml
name: i.hyper.smac Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Install GRASS GIS
      run: |
        sudo apt-get update
        sudo apt-get install grass grass-dev libradtran-dev
    - name: Run Tests
      run: |
        cd i.hyper.smac/testsuite
        python run_tests.py --verbose
```

### Test Coverage

```bash
# Generate coverage report
python -m coverage run run_tests.py
python -m coverage report
python -m coverage html
```

## Documentation

- **Programmer Manual**: `i_hyper_smac.dox` (Doxygen format)
- **API Reference**: Complete function documentation
- **Algorithm Details**: Mathematical formulations
- **Usage Examples**: Practical implementation examples

## License

Copyright (C) 2025 by the GRASS Development Team

This program is free software under the GNU General Public License (>=v2).
Read the file COPYING that comes with GRASS for details.

## Support

- **GRASS GIS Website**: https://grass.osgeo.org
- **Documentation**: https://grass.osgeo.org/documentation/
- **Mailing List**: grass-dev@lists.osgeo.org
- **Issue Tracker**: https://github.com/OSGeo/grass/issues
