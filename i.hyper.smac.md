# i.hyper.smac

## NAME

**i.hyper.smac** - Atmospheric correction of hyperspectral imagery using SMAC

## KEYWORDS

imagery, atmospheric correction, hyperspectral, reflectance, SMAC

## SYNOPSIS

**i.hyper.smac**

**i.hyper.smac --help**

**i.hyper.smac** input=*name* output=*name* dem=*name* solar_zenith=*float* [method=*string*] [sensor=*string*] [solar_azimuth=*float*] [view_zenith=*float*] [view_azimuth=*float*] [aod=*float*] [water_vapor=*float*] [ozone=*float*] [pressure=*float*] [aerosol_model=*string*] [visibility=*float*] [**-k**] [**-g**] [--overwrite] [--verbose] [--quiet]

## DESCRIPTION

**i.hyper.smac** performs atmospheric correction on hyperspectral imagery using the Simplified Method for Atmospheric Correction (SMAC). It converts top-of-atmosphere (TOA) radiance to surface reflectance by removing atmospheric effects.

The module corrects for:

- **Rayleigh scattering**: Molecular scattering that increases toward shorter wavelengths
- **Aerosol scattering and absorption**: Particle-related atmospheric effects
- **Gaseous absorption**: H2O, O3, O2, CO2, CH4, NO2, CO absorption bands

Two correction methods are available:

1. **simple**: Fast analytical correction using wavelength-dependent empirical formulas
2. **libradtran**: Accurate correction using pre-generated SMAC coefficients fitted to libRadtran radiative transfer simulations

### Algorithm

The surface reflectance (rho) is computed from TOA reflectance (rho*) using:

```
rho = (rho* - rho_atm) / (T_gas * T_scat * (1 + s*rho))
```

Where:
- rho_atm = atmospheric path reflectance
- T_gas = total gaseous transmission
- T_scat = scattering transmission (Rayleigh + aerosol)
- s = spherical albedo of the atmosphere

The SMAC coefficients encode these atmospheric functions as simple polynomial
formulas fitted to full radiative transfer simulations:

- **Gaseous transmission**: T = exp(a * (u*m)^n)
- **Total transmission**: T = a0T + a1T*tau/cos(theta) + (a2T*P + a3T)/(1+cos(theta))
- **Spherical albedo**: s = a0s*P + a3s + a1s*tau + a2s*tau^2
- **Rayleigh optical thickness**: Hansen & Travis (1974) analytical formula
- **Aerosol scaling**: Angstrom law tau(lambda) = tau(550) * (lambda/550)^(-alpha)

### Automatic Parameter Estimation

When atmospheric parameters are not provided, the module can estimate them:

- **AOD (Aerosol Optical Depth)**: Uses the Dense Dark Vegetation (DDV) method, identifying dark vegetation pixels and inferring AOD from blue band reflectance
- **Water Vapor**: Uses the differential absorption technique with bands near 940nm and reference bands at 865nm
- **Ozone**: Uses the Chappuis band method (500-700nm), converted from DU to cm-atm

## PARAMETERS

### Required Parameters

**input**=*name*
:   Name of input 3D raster map containing TOA radiance (W/m^2/sr/um)

**output**=*name*
:   Name for output 3D raster map (surface reflectance)

**dem**=*name*
:   Name of Digital Elevation Model raster (meters) for atmospheric pressure correction

**solar_zenith**=*float*
:   Solar zenith angle in degrees (0-90)

### Optional Parameters

**method**=*string*
:   Atmospheric correction method. Options: simple, libradtran. Default: simple

**sensor**=*string*
:   Sensor type for automatic band configuration. Options: PRISMA, AVIRIS, AVIRIS_NG, HYPERION, ENMAP, OSK_GHOST, PIXXEL, ESPER, IPERLITE, KUVASPACE_23, KUVASPACE_32, WYVERN_23, WYVERN_32, HYP4U, TANAGER

**solar_azimuth**=*float*
:   Solar azimuth angle in degrees (0-360). Default: 0

**view_zenith**=*float*
:   View (sensor) zenith angle in degrees (0-90). Default: 0 (nadir)

**view_azimuth**=*float*
:   View azimuth angle in degrees (0-360). Default: 0

**aod**=*float*
:   Aerosol optical depth at 550nm. If not specified, estimated automatically using DDV method

**water_vapor**=*float*
:   Water vapor column content in g/cm^2. If not specified, estimated from absorption bands

**ozone**=*float*
:   Ozone column content in cm-atm. Default: 0.3

**pressure**=*float*
:   Atmospheric pressure in hPa. If not specified, estimated from DEM

**aerosol_model**=*string*
:   Aerosol model type. Options: continental, maritime, urban, desert. Default: continental

**visibility**=*float*
:   Visibility in km. If not provided, estimated from AOD

### Flags

**-k**
:   Keep temporary bands (for debugging)

**-g**
:   Generate SMAC coefficients from libRadtran at runtime instead of using
    pre-generated files. Requires libRadtran and scipy. Only applies to
    method=libradtran. Generated coefficients are cached to the COEFS directory
    for reuse in subsequent runs.

## EXAMPLES

### Basic atmospheric correction with automatic parameter estimation

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 method=libradtran
```

### Correction with specified atmospheric parameters

```bash
i.hyper.smac input=prisma_toa output=prisma_sr dem=elevation \
    solar_zenith=28 solar_azimuth=145 method=libradtran \
    aod=0.12 water_vapor=1.8 ozone=0.32
```

### Maritime aerosol model for coastal scenes

```bash
i.hyper.smac input=coastal_toa output=coastal_sr dem=elevation \
    solar_zenith=40 method=libradtran \
    aerosol_model=maritime aod=0.08 water_vapor=3.5
```

### Off-nadir viewing geometry

```bash
i.hyper.smac input=offnadir_toa output=offnadir_sr dem=elevation \
    solar_zenith=30 solar_azimuth=120 \
    view_zenith=15 view_azimuth=270 \
    method=libradtran aod=0.2
```

### Runtime coefficient generation with -g flag

```bash
# Generate SMAC coefficients from libRadtran on the fly
# Useful when pre-generated files are not available or exact wavelength match is needed
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 method=libradtran aerosol_model=continental -g
```

## NOTES

### Input Data Requirements

- Input must be a 3D raster (GRASS RASTER3D) with spectral bands as the third dimension
- Values should be in TOA radiance (W/m^2/sr/um)
- Band wavelengths and FWHM should be stored in the raster metadata (as set by i.hyper.import)

### SMAC Coefficient Files

The libradtran method uses pre-generated SMAC coefficient files from the
`COEFS/` directory. Each file contains 49 coefficients (19 lines) encoding
the atmospheric correction model for one wavelength and aerosol type.

**Pre-generated coefficients**:
- Wavelength range: 350-2500nm at 10nm intervals (216 wavelengths)
- Aerosol models: Continental, Maritime, Urban, Desert
- Total: 864 coefficient files

At runtime, the module selects the nearest available coefficient file (max 5nm
mismatch with the 10nm grid). For exact wavelength matching, use the **-g** flag.

**Generating custom coefficients**:

The `smac_coef_generator.py` tool generates coefficients by running libRadtran
simulations and fitting the SMAC analytical formulas using `scipy.optimize.curve_fit`.

```bash
# Single wavelength
python3 lib/smac_coef_generator.py single 550 --aerosol continental --verbose

# Batch at 10nm resolution (default)
python3 lib/smac_coef_generator.py batch

# Batch at 1nm resolution for maximum spectral accuracy
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1

# Single aerosol type at 1nm (faster, ~2h)
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1 \
    --aerosol continental --verbose
```

### Aerosol Models

| Model | Description | Typical Conditions |
|-------|-------------|-------------------|
| Continental | Rural/background aerosol | Inland areas, moderate pollution |
| Maritime | Sea salt dominated | Coastal and oceanic areas |
| Urban | Soot and industrial particles | Cities, industrial regions |
| Desert | Mineral dust | Arid regions, dust events |

### Performance Considerations

- The simple method is faster but less accurate than the libradtran method
- The **-g** flag adds overhead (several seconds per band for libRadtran fitting) but generated coefficients are cached for reuse
- Automatic AOD/WVC estimation adds overhead but improves accuracy
- For time-series processing, consider pre-computing atmospheric parameters and pre-generating coefficients

### Limitations

- Assumes horizontally homogeneous atmosphere
- Best accuracy for AOD < 0.5 and SZA < 70 degrees
- May underperform in presence of clouds, shadows, or extreme aerosol loading
- Water bodies and snow may not be accurately corrected

## REFERENCES

Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143. https://doi.org/10.1080/01431169408954055

Hagolle, O. SMAC Python implementation. GitHub repository. https://github.com/olivierhagolle/SMAC

libRadtran radiative transfer package. https://www.libradtran.org/

Hansen, J.E. & Travis, L.D. (1974). Light scattering in planetary atmospheres. *Space Science Reviews*, 16, 527-610.

## SEE ALSO

*[i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html)* - Atmospheric correction using 6S model

*[i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html)* - Vegetation indices

*[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html)* - 3D raster information

## AUTHOR

Yann Chemin, dr.yann.chemin@gmail.com

*Based on SMAC algorithm by H. Rahman and G. Dedieu (CESBIO) and Python implementation by Olivier Hagolle*
