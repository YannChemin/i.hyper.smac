# i.hyper.smac - SMAC Atmospheric Correction for Hyperspectral Data

## DESCRIPTION
This module implements the Simplified Method for Atmospheric Correction (SMAC) for hyperspectral data in GRASS GIS. It provides two methods for atmospheric correction:

1. **Simple SMAC**: A faster, less accurate method suitable for quick processing
2. **Libradtran SMAC**: A more accurate method using the libRadtran radiative transfer model

The module is based on the original SMAC algorithm developed by [Rahman and Dedieu (1994)](https://doi.org/10.1016/0034-4257(94)90119-8) and the implementation by [Olivier Hagolle](https://github.com/olivierhagolle/SMAC).

## FEATURES
- Supports multiple hyperspectral sensors (AVIRIS, HYPERION, PRISMA, etc.)
- Two correction methods: simple and libradtran-based
- Automatic estimation of atmospheric parameters:
  - Aerosol Optical Depth (AOD) using the Dense Dark Vegetation (DDV) method
  - Water Vapor Content (WVC) using atmospheric absorption features
- Uses libRadtran for accurate atmospheric simulations (for libradtran method)
- Handles both per-pixel and scene-average atmospheric correction
- Preserves the original data structure and metadata

## REQUIREMENTS
- GRASS GIS 8.0 or later
- Python 3.x
- NumPy
- PyEphem (for solar zenith angle estimation)
- libRadtran (for atmospheric simulations)
- scipy (for interpolation and optimization)
- pandas (for data handling)
- scikit-learn (for regression models in AOD estimation)

## USAGE
```bash
# Simple SMAC (default)
i.hyper.smac input=name output=name elevation=name [sensor=string]
             [aod=float] [water_vapor=float] [ozone=float] 
             [solar_zenith=float] [solar_azimuth=float]
             [view_zenith=float] [view_azimuth=float]
             [--overwrite] [--verbose] [--quiet] [--keep-temp]

# Libradtran SMAC
i.hyper.smac input=name output=name elevation=name method=libradtran
             sensor=string [aod=float] [water_vapor=float] [ozone=float]
             [solar_zenith=float] [solar_azimuth=float]
             [view_zenith=float] [view_azimuth=float]
             [aot_lut=name] [visibility=float]
             [aerosol_type=string] [--overwrite] [--verbose] [--quiet] [--keep-temp]
```

## PARAMETERS
### Input/Output
- **input** - Name of input 3D raster map
- **output** - Name for output 3D raster map
- **elevation** - Name of elevation raster map (default: DEM from GRASS location)

### Method Selection
- **method** - Atmospheric correction method (simple, libradtran; default: simple)

### Atmospheric Parameters (both methods)
- **sensor** - Sensor type (AVIRIS, HYPERION, PRISMA, etc.)
- **aod** - Aerosol Optical Depth at 550nm (if not provided, will be estimated)
- **water_vapor** - Water vapor content in g/cm² (if not provided, will be estimated)
- **ozone** - Ozone content in cm-atm (default: 0.3)
- **solar_zenith** - Solar zenith angle in degrees (required)
- **solar_azimuth** - Solar azimuth angle in degrees (default: 0)
- **view_zenith** - View zenith angle in degrees (default: 0)
- **view_azimuth** - View azimuth angle in degrees (default: 0)
- **aerosol_type** - Aerosol type (continental, maritime, urban, desert; default: continental)

### Libradtran-specific Parameters
- **aot_lut** - Look-up table for AOD estimation (if not provided, will be generated)
- **visibility** - Visibility in km (if not provided, will be estimated from AOD)

## EXAMPLES
### Simple SMAC with automatic parameter estimation
```bash
# Automatically estimate AOD and water vapor
i.hyper.smac input=hyperspectral output=corrected elevation=dem sensor=AVIRIS
```

### Simple SMAC with manual parameter specification
```bash
# Manually specify AOD and water vapor
i.hyper.smac input=hyperspectral output=corrected elevation=dem sensor=AVIRIS \
  aod=0.2 water_vapor=2.5
```

### Simple SMAC with specific atmospheric parameters
```bash
i.hyper.smac input=hyperspectral output=corrected elevation=dem sensor=HYPERION \
  aod=0.2 water_vapor=2.5 solar_zenith=30
```

### Libradtran SMAC with custom parameters
```bash
i.hyper.smac input=hyperspectral output=corrected elevation=dem method=libradtran \
  sensor=AVIRIS aod=0.15 water_vapor=2.0 solar_zenith=25 \
  aerosol_type=maritime visibility=50
```

### Using an AOT look-up table with Libradtran
```bash
i.hyper.smac input=hyperspectral output=corrected elevation=dem method=libradtran \
  sensor=PRISMA aot_lut=aot_lookup.txt solar_zenith=30
```

## NOTES
- **Automatic Parameter Estimation**:
  - AOD is estimated using the Dense Dark Vegetation (DDV) method, which identifies dark vegetation pixels and uses their reflectance in the blue region
  - Water vapor is estimated using atmospheric absorption features around 940nm and 1130nm
  - If estimation fails, default values are used (AOD=0.15, Water Vapor=2.0 g/cm²)
  
- **Performance Considerations**:
  - The simple method is faster but less accurate than the libradtran method
  - AOD and WVC estimation add computational overhead but improve accuracy
  - The module creates temporary files in the current mapset (use `--keep-temp` to preserve them)

- **Requirements**:
  - The libradtran method requires the libRadtran software to be installed and properly configured
  - Additional Python packages (scipy, pandas, scikit-learn) are required for parameter estimation

- **Supported Sensors**:
  - PRISMA, AVIRIS, AVIRIS_NG, HYPERION, ENMAP, OSK_GHOST, PIXEL, ESPER, IPERLITE, KUVASPACE_23, KUVASPACE_32, WYVERN_23, WYVERN_32, HYP4U, TANAGER

- **Best Practices**:
  - For best results, provide as much metadata as possible
  - Verify the estimated parameters using the module's output messages
  - Consider using the `--keep-temp` flag for debugging if results are unexpected

## REFERENCES
- Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143.
- Original SMAC implementation: [https://github.com/olivierhagolle/SMAC](https://github.com/olivierhagolle/SMAC)

## AUTHORS
- Original SMAC algorithm: H. Rahman and G. Dedieu
- GRASS GIS implementation: [Yann Chemin] <dr.yann.chemin@gmail.com>

## SOURCE CODE
Available at: [GRASS GIS Addons Repository](https://github.com/OSGeo/grass-addons)

## SEE ALSO
- [r.smap](https://grass.osgeo.org/grass-stable/manuals/r.smap.html)
- [i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html)
- [i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html)
