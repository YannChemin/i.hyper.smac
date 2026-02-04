# i.hyper.smac

## NAME

**i.hyper.smac** - Atmospheric correction of hyperspectral imagery using SMAC

## KEYWORDS

imagery, atmospheric correction, hyperspectral, reflectance, SMAC

## SYNOPSIS

**i.hyper.smac**

**i.hyper.smac --help**

**i.hyper.smac** input=*name* output=*name* elevation=*name* solar_zenith=*float* [sensor=*string*] [solar_azimuth=*float*] [view_zenith=*float*] [view_azimuth=*float*] [aod=*float*] [water_vapor=*float*] [ozone=*float*] [aerosol_type=*string*] [coef_dir=*string*] [--overwrite] [--verbose] [--quiet]

## DESCRIPTION

**i.hyper.smac** performs atmospheric correction on hyperspectral imagery using the Simplified Method for Atmospheric Correction (SMAC). It converts top-of-atmosphere (TOA) radiance or reflectance to surface reflectance by removing atmospheric effects.

The module corrects for:

- **Rayleigh scattering**: Molecular scattering that increases toward shorter wavelengths (blue)
- **Aerosol scattering and absorption**: Particle-related atmospheric effects
- **Gaseous absorption**: H2O, O3, O2, CO2, CH4, NO2, CO absorption bands

The SMAC algorithm uses pre-computed coefficients derived from radiative transfer simulations. These coefficients encode the atmospheric correction parameters as simple polynomial functions, enabling fast processing of large hyperspectral datasets.

### Algorithm

The surface reflectance (ρ) is computed from TOA reflectance (ρ*) using:

```
ρ = (ρ* - ρ_atm) / (T_gas × T_scat × (1 + s×ρ))
```

Where:
- ρ_atm = atmospheric path reflectance
- T_gas = total gaseous transmission
- T_scat = scattering transmission (Rayleigh + aerosol)
- s = spherical albedo of the atmosphere

### Automatic Parameter Estimation

When atmospheric parameters are not provided, the module can estimate them:

- **AOD (Aerosol Optical Depth)**: Uses the Dense Dark Vegetation (DDV) method, identifying dark vegetation pixels and inferring AOD from blue band reflectance
- **Water Vapor**: Uses the differential absorption technique with bands near 940nm and reference bands at 865nm
- **Ozone**: Uses climatological values based on latitude and day of year, or can retrieve from TOMS/OMI data

## PARAMETERS

### Required Parameters

**input**=*name*
:   Name of input 3D raster map containing TOA radiance or reflectance

**output**=*name*
:   Name for output 3D raster map (surface reflectance)

**elevation**=*name*
:   Name of elevation raster map (DEM) for atmospheric pressure correction

**solar_zenith**=*float*
:   Solar zenith angle in degrees (0-90)

### Optional Parameters

**sensor**=*string*
:   Sensor type for automatic band configuration. Options: PRISMA, AVIRIS, AVIRIS_NG, HYPERION, ENMAP, EMIT, TANAGER, GENERIC_VNIR, GENERIC_FULL

**solar_azimuth**=*float*
:   Solar azimuth angle in degrees (0-360). Default: 0

**view_zenith**=*float*
:   View (sensor) zenith angle in degrees (0-90). Default: 0 (nadir)

**view_azimuth**=*float*
:   View azimuth angle in degrees (0-360). Default: 0

**aod**=*float*
:   Aerosol optical depth at 550nm. Range: 0.0-2.0. If not specified, estimated automatically using DDV method

**water_vapor**=*float*
:   Water vapor column content in g/cm². Range: 0.1-7.0. If not specified, estimated from absorption bands

**ozone**=*float*
:   Ozone column content in cm-atm. Range: 0.2-0.5. Default: 0.30

**aerosol_type**=*string*
:   Aerosol model type. Options: continental, maritime, urban, desert. Default: continental

**coef_dir**=*string*
:   Directory containing SMAC coefficient files. Default: module's COEFS directory

## EXAMPLES

### Basic atmospheric correction with automatic parameter estimation

```bash
i.hyper.smac input=aviris_toa output=aviris_sr elevation=dem \
    solar_zenith=35 sensor=AVIRIS
```

### Correction with specified atmospheric parameters

```bash
i.hyper.smac input=prisma_toa output=prisma_sr elevation=dem \
    solar_zenith=28 solar_azimuth=145 sensor=PRISMA \
    aod=0.12 water_vapor=1.8 ozone=0.32
```

### Maritime aerosol model for coastal scenes

```bash
i.hyper.smac input=coastal_toa output=coastal_sr elevation=dem \
    solar_zenith=40 sensor=AVIRIS_NG \
    aerosol_type=maritime aod=0.08 water_vapor=3.5
```

### Off-nadir viewing geometry

```bash
i.hyper.smac input=offnadir_toa output=offnadir_sr elevation=dem \
    solar_zenith=30 solar_azimuth=120 \
    view_zenith=15 view_azimuth=270 \
    sensor=HYPERION aod=0.2
```

### Using custom coefficient directory

```bash
i.hyper.smac input=custom_toa output=custom_sr elevation=dem \
    solar_zenith=25 coef_dir=/path/to/my/coefficients
```

## NOTES

### Input Data Requirements

- Input must be a 3D raster (GRASS RASTER3D) with spectral bands as the third dimension
- Values should be in TOA radiance (W/m²/sr/nm) or TOA reflectance (0-1)
- Band wavelengths should be stored in the raster metadata or specified via sensor type

### Coefficient Files

SMAC requires wavelength-specific coefficient files. The module includes pre-generated coefficients for:
- Wavelength range: 400-2500nm at 50nm intervals
- Aerosol types: Continental, Maritime, Urban, Desert
- Total: 172 coefficient files

For wavelengths not exactly matching available coefficients, the module uses the nearest available wavelength.

### Aerosol Models

| Model | Description | Typical Conditions |
|-------|-------------|-------------------|
| Continental | Rural/background aerosol | Inland areas, moderate pollution |
| Maritime | Sea salt dominated | Coastal and oceanic areas |
| Urban | Soot and industrial particles | Cities, industrial regions |
| Desert | Mineral dust | Arid regions, dust events |

### Performance Considerations

- The SMAC algorithm is computationally efficient, suitable for large hyperspectral datasets
- Automatic AOD/WVC estimation adds overhead but improves accuracy
- For time-series processing, consider pre-computing atmospheric parameters

### Limitations

- Assumes horizontally homogeneous atmosphere
- Best accuracy for AOD < 0.5 and SZA < 70°
- May underperform in presence of clouds, shadows, or extreme aerosol loading
- Water bodies and snow may not be accurately corrected

## REFERENCES

Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143. https://doi.org/10.1080/01431169408954055

Hagolle, O. SMAC Python implementation. GitHub repository. https://github.com/olivierhagolle/SMAC

## SEE ALSO

*[i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html)* - Atmospheric correction using 6S model

*[i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html)* - Vegetation indices

*[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html)* - 3D raster information

*[r3.stats](https://grass.osgeo.org/grass-stable/manuals/r3.stats.html)* - 3D raster statistics

## AUTHOR

Yann Chemin, dr.yann.chemin@gmail.com

*Based on SMAC algorithm by H. Rahman and G. Dedieu (CESBIO) and Python implementation by Olivier Hagolle*
