# i.hyper.smac - Apply SMAC atmospheric correction to hyperspectral imagery

## DESCRIPTION

`i.hyper.smac` applies atmospheric correction to hyperspectral imagery using the
SMAC (Simplified Method for Atmospheric Correction) algorithm enhanced with
libRadtran-based Look-Up Tables (LUTs). The module supports:

- **Full multiple-scattering** atmospheric correction via libRadtran DISORT
- **Parallel LUT generation** using 75% of available CPU resources
- **OpenCL acceleration** for atmospheric correction pixel processing
- **Spectral polishing** to remove outlier bands
- **Uncertainty quantification** for atmospheric correction results
- **Multiple sensor support** with pre-configured parameters

## KEY FEATURES

### Parallel LUT Generation (NEW)
- **75% Resource Usage**: Automatically uses 75% of available CPU resources
- **Smart Detection**: Auto-detects optimal worker count based on CPU cores and RAM
- **Memory Aware**: Considers system RAM constraints for uvspec processes
- **Performance**: 6x speedup on 8-core systems vs sequential
- **CPU-only**: uvspec runs on CPU with multiple parallel processes

### OpenCL Acceleration for Atmospheric Correction
- **Two-Stage Process**: LUT generation (CPU) + Atmospheric correction (OpenCL)
- **Device Selection**: Auto, GPU, or CPU OpenCL devices for correction stage
- **Memory Management**: Configurable memory limits for OpenCL kernels
- **Fallback**: Automatic CPU OpenCL fallback if GPU unavailable
- **Note**: uvspec itself is CPU-only; OpenCL accelerates pixel processing

## ARCHITECTURE

### Two-Stage Atmospheric Correction Process

**Stage 1: LUT Generation (CPU-only)**
```bash
# Multiple uvspec processes run in parallel on CPU
uvspec < input_1.inp  # Process 1: AOD=0.0, H2O=0.3, Albedo=0.0
uvspec < input_2.inp  # Process 2: AOD=0.1, H2O=0.3, Albedo=0.0
uvspec < input_3.inp  # Process 3: AOD=0.3, H2O=0.3, Albedo=0.0
# ... (96 total combinations)
```
- **libRadtran uvspec** calculates radiative transfer
- **Multiple CPU cores** run different parameter combinations
- **Generates complete lookup table** for all atmospheric conditions

**Stage 2: Atmospheric Correction (OpenCL)**
```bash
# OpenCL kernels apply LUT to hyperspectral data
opencl_kernel_process_pixels(lut_data, hyperspectral_image)
```
- **OpenCL kernels** process millions of pixels in parallel
- **GPU or CPU OpenCL** devices can be used
- **Applies pre-computed LUT** to actual scene data

## USAGE

### Basic Atmospheric Correction
```bash
i.hyper.smac input=hyperspectral output=corrected dem=elevation
```

### Advanced Configuration with Parallel Processing
```bash
i.hyper.smac input=hyperspectral output=corrected dem=elevation \
    solar_zenith=40 solar_azimuth=130 view_zenith=50 view_azimuth=100 \
    sensor=TANAGER aerosol_model=urban water_vapor=joint \
    parallel_lut=auto -c
```

### GPU Acceleration
```bash
i.hyper.smac input=hyperspectral output=corrected dem=elevation \
    opencl_device=gpu opencl_memory=2048
```

## PARAMETERS

### Required Parameters
- **input** - Input hyperspectral 3D raster map
- **output** - Output atmospherically corrected 3D raster map  
- **dem** - Digital Elevation Model in meters

### Atmospheric Parameters
- **aod** - Aerosol Optical Depth at 550nm (auto-estimated if not provided)
- **water_vapor** - Water vapor content (g/cm²) or estimation method
- **aerosol_model** - Aerosol type (continental, maritime, urban, desert)
- **pressure** - Atmospheric pressure (hPa)
- **ozone** - Ozone content (cm-atm, default: 0.3)

### Geometric Parameters
- **solar_zenith** - Solar zenith angle in degrees
- **solar_azimuth** - Solar azimuth angle in degrees (default: 0)
- **view_zenith** - View zenith angle in degrees (default: 0)
- **view_azimuth** - View azimuth angle in degrees (default: 0)

### Performance Parameters
- **parallel_lut** - Parallel LUT generation (auto, enabled, disabled)
  - **auto**: Use 75% of resources, auto-detect GPU/CPU (default)
  - **enabled**: Force parallel processing with 75% resources
  - **disabled**: Sequential processing (single core)
- **opencl_device** - OpenCL device type (auto, gpu, cpu)
- **opencl_memory** - OpenCL memory limit in MB (default: 1024)

### Sensor Configuration
- **sensor** - Pre-configured sensor parameters:
  - PRISMA, AVIRIS, AVIRIS_NG, HYPERION, ENMAP
  - OSK_GHOST, PIXXEL, ESPER, IPERLITE
  - KUVASPACE_23, KUVASPACE_32
  - WYVERN_23, WYVERN_32, HYP4U, TANAGER

### Advanced Options
- **adjacency_psf** - Adjacency effect PSF radius in km (default: 0)
- **output_uncertainty** - Uncertainty output raster (requires -u flag)

## FLAGS

- **-k** - Keep temporary bands
- **-p** - Apply spectral polishing to remove outlier bands
- **-c** - Clear cached LUT and regenerate from libRadtran
- **-u** - Compute per-band reflectance uncertainty
- **--overwrite** - Overwrite existing output files

## PERFORMANCE OPTIMIZATION

### Parallel LUT Generation
The module automatically uses 75% of available resources:

| System Cores | Workers Used | Speedup |
|--------------|--------------|---------|
| 4 cores      | 3 workers    | 3x      |
| 8 cores      | 6 workers    | 6x      |
| 16 cores     | 12 workers   | 12x     |

### Resource Usage Examples
```bash
# Default: 75% resources, auto-detect
i.hyper.smac input=rad output=smac dem=DEM parallel_lut=auto

# Force enable parallel processing
i.hyper.smac input=rad output=smac dem=DEM parallel_lut=enabled

# Disable for debugging
i.hyper.smac input=rad output=smac dem=DEM parallel_lut=disabled
```

### GPU Acceleration
```bash
# Auto-detect GPU
i.hyper.smac input=rad output=smac dem=DEM opencl_device=auto

# Force GPU usage
i.hyper.smac input=rad output=smac dem=DEM opencl_device=gpu

# Increase GPU memory
i.hyper.smac input=rad output=smac dem=DEM opencl_memory=2048
```

## EXAMPLES

### Basic Usage with Auto-Detection
```bash
# Automatic atmospheric correction with 75% parallel resources
i.hyper.smac input=PRISMA_L1C output=PRISMA_L2A dem=SRTM \
    sensor=PRISMA water_vapor=joint -c
```

### Advanced Configuration
```bash
# Full configuration with parallel processing and GPU acceleration
i.hyper.smac input=AVIRIS_NG output=AVIRIS_NG_corrected dem=DEM \
    solar_zenith=35.2 solar_azimuth=145.8 \
    view_zenith=15.5 view_azimuth=90.3 \
    aerosol_model=maritime water_vapor=1.2 \
    parallel_lut=enabled opencl_device=gpu \
    opencl_memory=2048 -p -u
```

### Batch Processing
```bash
# Process multiple scenes with optimal settings
for scene in scene1 scene2 scene3; do
    i.hyper.smac input=${scene}_raw output=${scene}_corrected dem=DEM \
        sensor=TANAGER parallel_lut=auto -c --overwrite
done
```

## NOTES

### LUT Generation
- **First Run**: Requires libRadtran LUT generation (~10-20 minutes with 75% parallel)
- **Subsequent Runs**: Uses cached LUTs (~1-2 minutes)
- **Force Regenerate**: Use -c flag to regenerate LUTs
- **Parallel Speed**: 75% resource usage provides optimal balance

### Memory Requirements
- **Sequential**: ~200MB per uvspec process
- **Parallel**: 75% of available RAM for LUT generation
- **GPU**: Additional 500MB-2GB for OpenCL processing

### Compatibility
- **GRASS GIS**: Version 8.0+
- **Python**: 3.8+
- **libRadtran**: 2.0+
- **OpenCL**: Optional for GPU acceleration

## AUTHORS

`i.hyper.smac` was developed by the GRASS GIS hyperspectral team with contributions
from the atmospheric correction and parallel computing communities.

## SEE ALSO

- `i.hyper.import` - Import hyperspectral data
- `i.hyper.spectral` - Spectral analysis tools
- `r.in.gdal` - Import raster data
- libRadtran documentation for atmospheric modeling

## COPYRIGHT

(C) 2024-2026 GRASS GIS Development Team

This program is free software under the GNU General Public License (>=v2).
Read the file COPYING that comes with GRASS for details.
