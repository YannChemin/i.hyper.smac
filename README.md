# i.hyper.smac - SMAC Atmospheric Correction for Hyperspectral Data

A GRASS GIS addon module implementing the Simplified Method for Atmospheric Correction (SMAC) for hyperspectral imagery.

## Overview

This module converts top-of-atmosphere (TOA) radiance/reflectance to surface reflectance by correcting for atmospheric effects including:
- Rayleigh (molecular) scattering
- Aerosol scattering and absorption
- Gaseous absorption (H2O, O3, O2, CO2, CH4, NO2, CO)

Based on the SMAC algorithm by Rahman & Dedieu (1994) and the implementation by Olivier Hagolle.

## Features

- **Hyperspectral Support**: Handles hundreds of spectral bands (400-2500nm)
- **Multiple Sensors**: Pre-configured for PRISMA, AVIRIS, AVIRIS-NG, HYPERION, EnMAP, EMIT, Tanager, and more
- **Automatic Parameter Estimation**:
  - Aerosol Optical Depth (AOD) using Dense Dark Vegetation (DDV) method
  - Water Vapor Content (WVC) using 940nm/1130nm absorption features
  - Ozone from climatological data or TOMS/OMI
- **Coefficient Generation**: Generate SMAC coefficients using libRadtran or analytical formulas
- **Four Aerosol Models**: Continental, Maritime, Urban, Desert

## Project Structure

```
i.hyper.smac/
├── i.hyper.smac.py          # Main GRASS GIS module
├── i.hyper.smac.html        # GRASS HTML manual
├── i.hyper.smac.md          # Module documentation
├── Makefile                 # GRASS module build file
├── lib/                     # Support library
│   ├── __init__.py          # Package exports
│   ├── smac.py              # Core SMAC algorithm
│   ├── aod.py               # AOD estimation (DDV method)
│   ├── wvc.py               # Water vapor estimation
│   ├── o3.py                # Ozone estimation
│   ├── radtran.py           # libRadtran integration
│   ├── smac_coef_generator.py  # Coefficient generation
│   └── utils.py             # Shared utilities
├── scripts/
│   └── generate_hyperspectral_coefs.py  # Batch coefficient generation
├── COEFS/                   # Pre-generated coefficient files (10nm, 350-2500nm)
│   ├── CONTINENTAL/         # Continental aerosol (216 files)
│   ├── MARITIME/            # Maritime aerosol (216 files)
│   ├── URBAN/               # Urban aerosol (216 files)
│   └── DESERT/              # Desert aerosol (216 files)
├── tests/                   # Unit tests
│   ├── test_smac.py         # SMAC algorithm tests
│   └── test_utils.py        # Utility function tests
└── docs/
    └── SMAC_COEFFICIENT_SPECIFICATION.md
```

## Requirements

### Required
- GRASS GIS 8.0+
- Python 3.8+
- NumPy
- SciPy

### Optional
- libRadtran (for coefficient generation from radiative transfer simulations)
- ephem (for solar position calculations)
- scikit-learn (for advanced AOD estimation)

## Installation

### As GRASS GIS Addon
```bash
# From GRASS GIS
g.extension extension=i.hyper.smac url=/path/to/i.hyper.smac
```

### Manual Installation
```bash
# Clone/copy to GRASS addons directory
cp -r i.hyper.smac $GISBASE/scripts/

# Or add to GRASS_ADDON_PATH
export GRASS_ADDON_PATH=/path/to/i.hyper.smac:$GRASS_ADDON_PATH
```

## Usage

### Basic Usage
```bash
# With automatic parameter estimation
i.hyper.smac input=toa_radiance output=surface_reflectance \
    elevation=dem solar_zenith=30 sensor=AVIRIS

# With manual atmospheric parameters
i.hyper.smac input=toa_radiance output=surface_reflectance \
    elevation=dem solar_zenith=30 sensor=PRISMA \
    aod=0.15 water_vapor=2.0 ozone=0.34
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| input | Input 3D raster (TOA radiance/reflectance) | required |
| output | Output 3D raster (surface reflectance) | required |
| elevation | DEM raster for pressure correction | required |
| sensor | Sensor type for band configuration | - |
| solar_zenith | Solar zenith angle (degrees) | required |
| solar_azimuth | Solar azimuth angle (degrees) | 0 |
| view_zenith | View zenith angle (degrees) | 0 |
| view_azimuth | View azimuth angle (degrees) | 0 |
| aod | Aerosol optical depth at 550nm | auto |
| water_vapor | Water vapor column (g/cm²) | auto |
| ozone | Ozone column (cm-atm) | 0.30 |
| aerosol_type | continental/maritime/urban/desert | continental |

### Supported Sensors

- PRISMA (ASI)
- AVIRIS, AVIRIS-NG (NASA)
- HYPERION (EO-1)
- EnMAP (DLR)
- EMIT (NASA/JPL)
- Tanager (Planet)
- Generic VNIR/Full-range configurations

## SMAC Coefficient Generation

SMAC uses 49 coefficients per spectral band to parameterize atmospheric
transmission, scattering, and absorption. These coefficients are generated
by fitting the SMAC analytical formulas to libRadtran radiative transfer
simulations across a range of atmospheric conditions (varying AOD, SZA,
surface albedo, gas amounts).

The fitted parameters include:
- **Gaseous absorption** (H2O, O3): `T = exp(a * (u*m)^n)`
- **Transmission**: `T = a0T + a1T*tau/cos(theta) + (a2T*P + a3T)/(1+cos(theta))`
- **Spherical albedo**: `s = a0s*P + a3s + a1s*tau + a2s*tau^2`
- **Rayleigh optical thickness** (analytical, Hansen & Travis 1974)
- **Aerosol scaling** (Angstrom law): `tau(lambda) = tau(550) * (lambda/550)^(-alpha)`

### Pre-generated Coefficients

The `COEFS/` directory ships with coefficients at 10nm intervals from 350-2500nm
for all four aerosol models (864 files total):

```
COEFS/
├── CONTINENTAL/   coef_350nm_CONTINENTAL.dat ... coef_2500nm_CONTINENTAL.dat
├── MARITIME/      coef_350nm_MARITIME.dat    ... coef_2500nm_MARITIME.dat
├── URBAN/         coef_350nm_URBAN.dat       ... coef_2500nm_URBAN.dat
└── DESERT/        coef_350nm_DESERT.dat      ... coef_2500nm_DESERT.dat
```

At runtime, `find_coef_file()` selects the nearest available wavelength file
(max 5nm mismatch with the 10nm grid).

### Generating Coefficients with `smac_coef_generator.py`

The generator runs libRadtran simulations and uses `scipy.optimize.curve_fit`
to fit the SMAC formulas. It requires **libRadtran** and **scipy**.

#### Single wavelength

```bash
# Generate coefficients for 550nm, continental aerosol
python3 lib/smac_coef_generator.py single 550 --aerosol continental --verbose

# Specify output file and FWHM
python3 lib/smac_coef_generator.py single 940 --fwhm 5 --output coef_940nm.dat --verbose
```

#### Batch generation (10nm resolution)

```bash
# Default: 350-2500nm at 10nm steps, all four aerosol types
python3 lib/smac_coef_generator.py batch

# Custom range and step
python3 lib/smac_coef_generator.py batch --min 400 --max 1000 --step 10

# Single aerosol type
python3 lib/smac_coef_generator.py batch --aerosol continental --verbose

# Custom output directory
python3 lib/smac_coef_generator.py batch --output-dir /path/to/coefs
```

#### Batch generation at 1nm resolution

For maximum spectral accuracy (e.g., narrow-band sensors or resolving
sharp gas absorption features), generate at 1nm intervals. This produces
2151 files per aerosol type (8604 total) and takes several hours:

```bash
# All aerosol types, 1nm steps (350-2500nm)
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1

# Single aerosol type for faster generation (~2h instead of ~8h)
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1 \
    --aerosol continental --verbose

# Only the VNIR range at 1nm
python3 lib/smac_coef_generator.py batch --min 400 --max 1000 --step 1

# Run in background with progress logging
nohup python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1 \
    > coef_generation.log 2>&1 &
tail -f coef_generation.log
```

### Runtime Generation with the `-g` Flag

Instead of using pre-generated files, coefficients can be generated on the
fly from libRadtran during atmospheric correction. Generated coefficients
are cached to the `COEFS/` directory for reuse:

```bash
i.hyper.smac input=toa_radiance output=surface_reflectance \
    dem=elevation method=libradtran solar_zenith=30 \
    aerosol_model=continental -g
```

This is useful when:
- Working at wavelengths not covered by pre-generated files
- Needing exact-wavelength coefficients without spectral interpolation
- Testing different aerosol models without pre-generating all files

### Python API

```python
from lib.smac_coef_generator import generate_smac_coefficients, generate_batch

# Generate a single coefficient set
coef = generate_smac_coefficients(
    wavelength_nm=550.0,
    fwhm_nm=10.0,
    aerosol_type='continental',
    output_file='coef_550nm.dat',
    verbose=True
)

# Batch generation
generate_batch(
    wavelength_min=400,
    wavelength_max=2500,
    step=10,
    aerosol_types=['continental', 'maritime'],
    output_dir='COEFS/'
)
```

## Library API

The `lib/` package can be used independently:

```python
from lib import smac, estimate_aod, estimate_wvc, estimate_ozone
from lib.smac_coef_generator import generate_smac_coefficients

# Load coefficients and apply SMAC
coef = smac.coeff('COEFS/CONTINENTAL/coef_550nm_CONTINENTAL.dat')
surface_refl = smac.smac_inv(toa_refl, theta_s, phi_s, theta_v, phi_v,
                              pressure, taup550, uo3, uh2o, coef)

# Estimate atmospheric parameters
aod = estimate_aod(red_band, nir_band, blue_band)
wvc = estimate_wvc(band_940nm, band_865nm)
o3 = estimate_ozone(latitude, day_of_year)

# Generate coefficients
coef = generate_smac_coefficients(
    wavelength_nm=550.0,
    fwhm_nm=10.0,
    aerosol_type='continental',
    output_file='coef_550nm.dat'
)
```

## Testing

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test module
python -m pytest tests/test_smac.py -v
```

Current test coverage: 23 tests covering:
- SMAC forward/inverse transforms
- Pressure altitude correction
- Round-trip accuracy
- Utility functions
- Atmospheric parameter estimation bounds

## Algorithm Details

### SMAC Equation

Surface reflectance ρ is computed from TOA reflectance ρ* using:

```
ρ = (ρ* - ρ_atm) / (T_gas × T_scat × (1 + s×ρ))
```

Where:
- ρ_atm = atmospheric path reflectance
- T_gas = gaseous transmission (H2O, O3, etc.)
- T_scat = scattering transmission (Rayleigh + aerosol)
- s = spherical albedo of atmosphere

### Gaseous Transmission
```
T_gas = exp(a × (u × m)^n)
```
- u = gas column content
- m = air mass (1/cos(θs) + 1/cos(θv))
- a, n = fitted coefficients

### Rayleigh Optical Thickness
```
τ_r = 0.008569 × λ^(-4) × (1 + 0.0113×λ^(-2))
```
Where λ is wavelength in micrometers.

## References

1. Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143.

2. Hagolle, O. SMAC Python implementation. https://github.com/olivierhagolle/SMAC

3. libRadtran radiative transfer package. https://www.libradtran.org/

## Authors

- **Original SMAC Algorithm**: H. Rahman and G. Dedieu (CESBIO)
- **GRASS GIS Implementation**: Yann Chemin <dr.yann.chemin@gmail.com>

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org/>

## See Also

- [i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html) - 6S-based atmospheric correction
- [i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html) - Vegetation indices
- [r3.in.xyz](https://grass.osgeo.org/grass-stable/manuals/r3.in.xyz.html) - Import 3D raster data
