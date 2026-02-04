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
├── COEFS/                   # Pre-generated coefficient files
│   ├── CONTINENTAL/         # Continental aerosol (43 files)
│   ├── MARITIME/            # Maritime aerosol (43 files)
│   ├── URBAN/               # Urban aerosol (43 files)
│   └── DESERT/              # Desert aerosol (43 files)
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

## Coefficient Generation

### Using Pre-generated Coefficients
The `COEFS/` directory contains coefficients for 400-2500nm at 50nm intervals for all four aerosol types (172 files total).

### Generate New Coefficients

#### With libRadtran (recommended)
```bash
cd scripts
python generate_hyperspectral_coefs.py --sensor PRISMA --aerosol continental
```

#### Analytical Mode (no libRadtran required)
```bash
python generate_hyperspectral_coefs.py --sensor AVIRIS --analytical
```

#### Custom Wavelength Range
```bash
python generate_hyperspectral_coefs.py --start 400 --end 1000 --step 5 --fwhm 5
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

GNU General Public License v2 or later (GPLv2+)

## See Also

- [i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html) - 6S-based atmospheric correction
- [i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html) - Vegetation indices
- [r3.in.xyz](https://grass.osgeo.org/grass-stable/manuals/r3.in.xyz.html) - Import 3D raster data
