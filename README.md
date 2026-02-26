# i.hyper.smac - SMAC Atmospheric Correction for Hyperspectral Data

A GRASS GIS addon module implementing atmospheric correction for hyperspectral imagery, combining the Simplified Method for Atmospheric Correction (SMAC) with a libRadtran look-up table (LUT) approach and ISOFIT-inspired enhancements.

**Current Implementation**: The module currently uses the LUT method (libRadtran DISORT) as the default and only implemented method. The simple and libradtran methods are documented for reference but are not yet implemented in this version.

## Overview

This module converts top-of-atmosphere (TOA) radiance to surface reflectance by correcting for atmospheric effects including:
- Rayleigh (molecular) scattering
- Aerosol scattering and absorption
- Gaseous absorption (H2O, O3, O2, CO2, CH4, NO2, CO)
- Adjacency effects (photon scattering from neighbouring pixels)

Based on the SMAC algorithm by Rahman & Dedieu (1994) and the implementation by Olivier Hagolle, extended with a libRadtran LUT inversion and techniques from the ISOFIT optimal estimation framework (Thompson et al., 2018).

## Features

- **LUT Method** (current implementation): Full multiple-scattering LUT from libRadtran DISORT with per-pixel H2O retrieval and ISOFIT-inspired enhancements
- **Simple Method** (planned): Fast analytical correction using empirical formulas  
- **libradtran Method** (planned): SMAC coefficients fitted to libRadtran simulations
- **Hyperspectral Support**: Handles hundreds of spectral bands (400-2500nm)
- **Multiple Sensors**: Pre-configured for PRISMA, AVIRIS, AVIRIS-NG, HYPERION, EnMAP, Tanager, and more
- **OpenCL GPU Acceleration**: GPU and multi-core acceleration for faster processing
- **Automatic Parameter Estimation**:
  - Aerosol Optical Depth (AOD) using MODIS Dark Target algorithm (Kaufman 1997, Levy 2007) with Dark Pixel and constant fallbacks
  - Water Vapor Content (WVC) using continuum-removal band-depth at 940nm/1130nm with air mass correction, multi-band averaging, and weighted combination
  - Ozone from: Chappuis band method (500-700nm)
- **ISOFIT-Inspired Enhancements** (LUT method):
  - Superpixel atmospheric retrieval with spatial interpolation
  - In-loop adjacency effect correction (Vermote et al., 1997)
  - Surface reflectance prior (3-component Gaussian mixture)
  - Per-band reflectance uncertainty propagation
  - MAP inner loop (prior-weighted inversion)
  - Model discrepancy matrix (per-band RT error floor)
- **Aerosol Spectral Calibration**:
  - Native Shettle/Fenn aerosol spectral shape (no Angstrom override)
  - Scene-retrieved Angstrom exponent from Dark Target
  - NIR transmittance correction from dark vegetation targets
- **SWIR H2O Accuracy**: Dense 7-point H2O LUT grid with log-space transmittance interpolation for accurate continuum absorption
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
│   ├── lut.py               # Atmospheric LUT (libRadtran DISORT)
│   ├── adjacency.py         # Adjacency effect correction (Vermote 1997)
│   ├── superpixel.py        # SLIC segmentation and spatial interpolation
│   ├── surface_model.py     # Surface prior (Gaussian mixture, MAP)
│   ├── uncertainty.py       # Reflectance uncertainty and model discrepancy
│   ├── spectral_polish.py   # Spectral outlier detection and smoothing
│   ├── aod.py               # AOD estimation (Dark Target method)
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
│   ├── DESERT/              # Desert aerosol (216 files)
│   └── LUT/                 # Cached atmospheric LUTs (.npz)
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
- libRadtran (for LUT generation and coefficient fitting)
- scikit-image (for superpixel segmentation — SLIC)
- ephem (for solar position calculations)

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

Atmospheric parameters (AOD, water vapour, ozone) are auto-estimated from
the scene when not specified. The aerosol model selects the Shettle/Fenn
spectral shape used for LUT generation; the atmosphere profile affects the
gas transmittance columns.

### Desert / arid land (Sahara, Arabian Peninsula)

Mineral dust has a coarse size distribution (low Ångström exponent) and
very low precipitable water. A wide AOD range is needed to cover heavy
dust events.

```bash
i.hyper.smac input=algeria_radiance output=algeria_refl \
    dem=elevation solar_zenith=40 solar_azimuth=160 \
    view_zenith=5 view_azimuth=310 \
    aerosol_model=desert \
    aod=0.5 water_vapor=0.4 ozone=0.28
```

### Tropical rainforest (Amazon, Congo)

Very high water vapour (> 5 g/cm²) and low AOD. Auto water vapour
retrieval from the 940 nm band works well in clean-sky tropical
conditions. Adjacency correction is important at forest-river boundaries.

```bash
i.hyper.smac input=amazon_radiance output=amazon_refl \
    dem=elevation solar_zenith=28 solar_azimuth=90 \
    view_zenith=0 view_azimuth=0 \
    aerosol_model=continental \
    adjacency_psf=1.5 -p
```

### Continental urban temperate mid-winter (e.g. Paris, January)

High SZA (> 70°), low H₂O, and elevated urban AOD. The `-u` flag is
particularly valuable here because the long atmospheric path inflates
AOD sensitivity. Auto-estimated AOD from the Dark Target algorithm
detects pollution over urban surfaces.

```bash
i.hyper.smac input=paris_radiance output=paris_refl \
    dem=elevation solar_zenith=72 solar_azimuth=175 \
    view_zenith=5 view_azimuth=55 \
    aerosol_model=urban \
    water_vapor=0.5 ozone=0.34 \
    output_uncertainty=paris_refl_unc -u
```

### Coastal temperate mid-summer (e.g. Mediterranean, July)

Low AOD, maritime sea-salt aerosol, moderate H₂O. Adjacency correction
is critical at land-sea boundaries: bright land reflectance leaks into
adjacent water pixels through the diffuse transmittance component of
maritime aerosol.

```bash
i.hyper.smac input=med_radiance output=med_refl \
    dem=elevation solar_zenith=22 solar_azimuth=200 \
    view_zenith=0 view_azimuth=0 \
    aerosol_model=maritime \
    adjacency_psf=0.5 \
    output_uncertainty=med_refl_unc -u
```

### Sub-arctic mid-winter — boreal forest under snow

Very low H₂O, low AOD, extremely high SZA. The MAP surface prior
(`-u` + surface regularisation in the LUT pipeline) constrains
bright-snow spectral spikes from residual gas-correction errors.

```bash
i.hyper.smac input=sweden_radiance output=sweden_refl \
    dem=elevation solar_zenith=78 solar_azimuth=175 \
    view_zenith=8 view_azimuth=355 \
    aerosol_model=continental \
    water_vapor=0.3 ozone=0.34 \
    output_uncertainty=sweden_refl_unc -u -p
```

### Sub-arctic mid-summer — tundra / boreal

Moderate sun, clean sub-arctic air with biogenic terpene aerosol.
Auto-estimated parameters work well in this low-AOD regime.

```bash
i.hyper.smac input=yukon_radiance output=yukon_refl \
    dem=elevation solar_zenith=45 solar_azimuth=180 \
    view_zenith=2 view_azimuth=0 \
    aerosol_model=continental -p
```

### Full ISOFIT-Inspired Processing

Combines all six enhancements: per-pixel AOD/H₂O maps, Gaussian spatial
smoothing, adjacency correction, surface prior MAP, uncertainty output,
and model discrepancy floor.

```bash
i.hyper.smac input=toa_radiance output=surface_reflectance \
    dem=elevation solar_zenith=30 solar_azimuth=180 \
    view_zenith=0 view_azimuth=0 adjacency_psf=1.0 \
    output_uncertainty=reflectance_unc -p -u
```

### Legacy Methods

```bash
# libradtran method with runtime coefficient generation
i.hyper.smac input=toa_radiance output=surface_reflectance \
    dem=elevation solar_zenith=30 method=libradtran -g

# Simple analytical method
i.hyper.smac input=toa_radiance output=surface_reflectance \
    dem=elevation solar_zenith=30 method=simple
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| input | Input 3D raster (TOA radiance W/m²/sr/um) | required |
| output | Output 3D raster (surface reflectance) | required |
| dem | DEM raster for pressure correction | required |
| solar_zenith | Solar zenith angle (degrees) | required |
| solar_azimuth | Solar azimuth angle (degrees) | required |
| view_zenith | View zenith angle (degrees) | required |
| view_azimuth | View azimuth angle (degrees) | required |
| method | Correction method: lut (only implemented) | lut |
| sensor | Sensor type for band configuration | - |
| aod | Aerosol optical depth at 550nm | auto |
| water_vapor | Water vapor column (g/cm²) | auto |
| ozone | Ozone column (cm-atm) | 0.30 |
| aerosol_model | continental/maritime/urban/desert | continental |
| angstrom | Angstrom exponent override (default: scene-retrieved or native Shettle) | auto |
| adjacency_psf | Adjacency PSF radius in km (0=disabled) | 0 |
| output_uncertainty | Output uncertainty 3D raster (requires -u) | - |
| opencl_device | OpenCL device type (auto,gpu,cpu) | auto |
| opencl_memory | OpenCL memory limit in MB (0=unlimited) | 1024 |

### Flags

| Flag | Description |
|------|-------------|
| -k | Keep temporary bands (debugging) |
| -g | Generate SMAC coefficients from libRadtran at runtime (planned feature) |
| -p | Apply spectral polishing (outlier band removal) |
| -c | Clear cached LUT and regenerate |
| -u | Compute per-band reflectance uncertainty |

### Supported Sensors

- PRISMA (ASI)
- AVIRIS, AVIRIS-NG (NASA)
- HYPERION (EO-1)
- EnMAP (DLR)
- Tanager (Planet)
- OSK GHOST, PIXXEL, ESPER, IPERLITE
- KUVASPACE (23/32 bands), WYVERN (23/32 bands), HYP4U

## OpenCL GPU Acceleration

The module supports OpenCL GPU and multi-core acceleration for significantly faster processing:

### Performance Benefits

- **3-15x speedup** depending on hardware and dataset size
- **Parallel processing**: Multiple bands processed simultaneously
- **Memory optimization**: Configurable memory limits for different hardware
- **Automatic fallback**: Graceful degradation to CPU processing

### Usage Examples

```bash
# Auto-detect best device (GPU preferred)
i.hyper.smac input=rad output=corr opencl_device=auto ...

# Force GPU acceleration
i.hyper.smac input=rad output=corr opencl_device=gpu ...

# Use CPU OpenCL (multi-core)
i.hyper.smac input=rad output=corr opencl_device=cpu ...

# Limit memory usage
i.hyper.smac input=rad output=corr opencl_memory=2048 ...
```

### Installation

```bash
# Install pyOpenCL
sudo apt install python3-pyopencl

# Verify installation
python3 -c "import pyopencl; print('OpenCL available')"
```

### Supported Devices

- **NVIDIA GPUs**: RTX series, GTX series, Tesla series
- **AMD GPUs**: RX series, Radeon series
- **Intel GPUs**: Arc series, integrated graphics
- **Multi-core CPUs**: Intel, AMD processors

See [OPENCL_INTEGRATION.md](OPENCL_INTEGRATION.md) for detailed documentation.

## LUT Atmospheric Correction

The LUT method (default) uses a precomputed 3D look-up table from libRadtran DISORT 16-stream simulations. The LUT is indexed by (wavelength, AOD, H2O) and stores four atmospheric parameters: path reflectance (R_atm), downward transmittance (T_down), upward transmittance (T_up), and spherical albedo (s).

**Alternative backend**: The companion module [i.hyper.atcorr](../i.hyper.atcorr/) implements the same LUT approach using a native C port of 6SV2.1 (no libRadtran required). It supports all six standard atmosphere profiles, four aerosol models, custom Mie aerosols, and nine BRDF surface models. The Python bindings (`lib/atcorr.py` + `lib/libatcorr.so`) allow i.hyper.smac to use the 6SV2.1 engine without a separate GRASS invocation.

### LUT Generation

Each LUT is built from two-albedo runs (albedo=0 and albedo=0.5) across 9 AOD values and 7 H2O values = 126 libRadtran simulations. The dense 7-point H2O grid (0.3x-2.5x scene mean) accurately captures the nonlinear H2O continuum absorption in SWIR bands (1500-2500nm). Transmittance parameters (T_down, T_up) are interpolated in log-space, since ln(T) is nearly linear in H2O column, eliminating the systematic bias from linear interpolation on the convex exp(-tau) function. The 6S-style decomposition separates direct and diffuse transmittance components.

Generated LUTs are cached as `.npz` files in `COEFS/LUT/` and reused when scene geometry matches. Use `-c` to force regeneration.

### Per-Pixel H2O Retrieval

Water vapor is retrieved per-pixel from the 940nm absorption band (ISOFIT-style). For each H2O value in the LUT grid, the observed radiance at 865nm (shoulder), 945nm (absorption), and 1040nm (shoulder) is inverted. The correct H2O produces a smooth continuum across the absorption feature.

### ISOFIT-Inspired Processing Pipeline

When enabled, the LUT method applies the following processing chain:

```
1. Load per-pixel maps (AOD, WVC, ozone, DEM)
2. Generate/load 3D LUT (7-point H2O grid)
3. SLIC superpixel segmentation on NIR band          [#1 superpixel]
4. Smooth AOD/H2O via superpixels                    [#1 superpixel]
5. Calibrate path radiance + NIR transmittance       [dark target calibration]
6. Retrieve per-pixel H2O from 940nm
7. Per-band loop:
   a. Extract band, convert to TOA reflectance
   b. LUT interpolation with per-pixel H2O           [log-space T interp]
   c. Path radiance calibration (VIS, <1000nm)
   d. NIR/SWIR transmittance correction (>=700nm)
   e. Algebraic inversion
   f. In-loop adjacency correction                   [#2 adjacency]
   g. Per-band uncertainty computation                [#4 uncertainty]
8. After all bands:
   a. Surface prior regularization / MAP blend        [#3/#5 surface prior]
   b. Model discrepancy added to uncertainty          [#6 discrepancy]
   c. Spectral polishing (on regularized data)
   d. Write reflectance + uncertainty rasters
```

## SMAC Coefficient Generation

SMAC uses 49 coefficients per spectral band to parameterize atmospheric
transmission, scattering, and absorption. These coefficients are generated
by fitting the SMAC analytical formulas to libRadtran radiative transfer
simulations across a range of atmospheric conditions.

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
```

#### Batch generation at 1nm resolution

For maximum spectral accuracy (e.g., narrow-band sensors or resolving
sharp gas absorption features), generate at 1nm intervals:

```bash
# All aerosol types, 1nm steps (350-2500nm)
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1

# Single aerosol type for faster generation (~2h instead of ~8h)
python3 lib/smac_coef_generator.py batch --min 350 --max 2500 --step 1 \
    --aerosol continental --verbose
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

## Library API

The `lib/` package can be used independently:

```python
from lib import smac, estimate_aod, estimate_wvc
from lib.lut import AtmosphericLUT
from lib.surface_model import SurfaceModel
from lib.uncertainty import compute_reflectance_uncertainty, compute_model_discrepancy

# LUT-based correction
lut = AtmosphericLUT.get_or_generate(
    sza=30, vza=0, phi=0, pressure=1013.25,
    aerosol_model='continental', h2o=1.5, o3=0.3
)
R_atm, T_down, T_up, s = lut.interpolate(wavelength=550, aod=0.15, h2o=1.5)

# Surface prior regularization
model = SurfaceModel(wavelengths)
regularized = model.regularize(reflectance_stack, weight=0.1)

# Model discrepancy
sigma_model = compute_model_discrepancy(wavelengths)

# Load coefficients and apply SMAC
coef = smac.coeff('COEFS/CONTINENTAL/coef_550nm_CONTINENTAL.dat')
surface_refl = smac.smac_inv(toa_refl, theta_s, phi_s, theta_v, phi_v,
                              pressure, taup550, uo3, uh2o, coef)

# Estimate AOD (Dark Target with fallback chain)
aod_map, aod = estimate_aod(
    input_raster, dem, method='auto',
    solar_zenith=30.0, view_zenith=0.0,
    solar_azimuth=180.0, view_azimuth=0.0,
    pressure=1013.25
)

# Estimate WVC (air mass-corrected, multi-band, weighted)
wvc_map, wvc = estimate_wvc(
    input_raster, dem, method='average',
    solar_zenith=30.0, view_zenith=0.0
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

Surface reflectance rho is computed from TOA reflectance rho* using:

```
rho = (rho* - rho_atm) / (T_gas * T_scat * (1 + s*rho))
```

Where:
- rho_atm = atmospheric path reflectance
- T_gas = gaseous transmission (H2O, O3, etc.)
- T_scat = scattering transmission (Rayleigh + aerosol)
- s = spherical albedo of atmosphere

### LUT Inversion (6S-style)

The LUT method decomposes the atmosphere using two-albedo simulations:

```
R_atm = pi * L_toa(albedo=0) / E0     # path reflectance
T_down = (E_dir + E_dif) / E0_toa      # downward transmittance
s = (E_global(rho1) / E_global(0) - 1) / (rho1 * E_global(rho1) / E_global(0))
T_up = pi * (L(rho1) - L(0)) * (1 - s*rho1) / (E0 * rho1 * T_down)
```

The inversion solves:
```
y = (rho_toa - R_atm) / (T_down * T_up)
rho_boa = y / (1 + s * y)
```

### Adjacency Effect Correction

Following Vermote et al. (1997), the total transmittance is decomposed into
direct (Beer-Lambert) and diffuse components. The diffuse component picks up
signal from the environmental reflectance (spatial average) rather than the
pixel reflectance:

```
T_diff = clip(T_down*T_up - T_dir, 0)
r_env = spatial_average(r_boa, PSF_radius)
correction = T_diff * s * (r_boa - r_env) / (1 - s*r_env)
r_adj = r_boa + correction
```

### Aerosol Spectral Calibration

Three techniques reduce systematic NIR/SWIR reflectance bias:

1. **Native Shettle/Fenn spectral shape**: Aerosol models use their internally
   consistent spectral properties from libRadtran's `aerosol_haze` directive,
   preserving the physically-based spectral variation of AOD.

2. **Scene-retrieved Angstrom exponent**: The Dark Target algorithm derives the
   Angstrom exponent from the blue/red AOD ratio. This is used as the LUT
   parameter when no explicit override is provided.

3. **NIR transmittance correction**: Dark vegetation targets at 860nm (known
   surface reflectance ~ 3.5 x rho_SWIR) provide an empirical transmittance
   correction applied to all bands >= 700nm:
```
t_corr = rho_retrieved(860nm) / rho_expected(860nm)
T_eff = T_down * T_up * t_corr
```

### SWIR H2O Continuum Accuracy

The 7-point H2O grid with log-space transmittance interpolation addresses the
nonlinear absorption in SWIR bands. Since T = exp(-tau) and tau is roughly
proportional to H2O column, interpolating ln(T) rather than T eliminates the
convexity bias that linear interpolation introduces. Testing shows this reduces
transmittance interpolation error from up to 54% (old 3-point linear) to <0.01%
at strong SWIR absorption bands.

### Surface Prior and MAP Estimate

The surface prior is a 3-component Gaussian mixture (vegetation, soil, water)
with per-band diagonal covariance. The MAP estimate blends the observation with
the prior:

```
r_map = (r_obs/sigma_obs^2 + r_prior/sigma_prior^2) / (1/sigma_obs^2 + 1/sigma_prior^2)
```

The measurement variance sigma_obs^2 includes instrument noise, AOD uncertainty,
and model discrepancy. In bands with high uncertainty (e.g., gas absorption
edges), the prior dominates; in clean bands, the observation dominates.

### Uncertainty Propagation

Per-band reflectance uncertainty combines three sources in quadrature:

1. **Instrument noise**: NEDL estimated from dark pixels, propagated as
   sigma_rfl = pi * NEDL * d^2 / (E0 * cos(SZA) * T_down * T_up)

2. **AOD uncertainty**: +/-0.04 AOD perturbation through the LUT inversion

3. **Model discrepancy**: Per-band systematic RT error floor, empirically
   characterized (0.5% baseline, up to 2.5% at gas band edges)

### AOD Estimation (Dark Target Algorithm)

Implements the MODIS Dark Target algorithm (Kaufman 1997, Levy 2007):

1. **Band selection**: Blue (470nm), Red (660nm), NIR (860nm), SWIR (2130nm)
2. **Dark Target pixel selection**: NDVI > 0.1, SWIR in 0.01-0.25, 20th-50th percentile
3. **Surface reflectance prediction**: rho_blue = 0.25 x rho_SWIR, rho_red = 0.50 x rho_SWIR
4. **Single-scattering AOD inversion** with Henyey-Greenstein phase function
5. **Angstrom exponent** from blue/red AOD ratio, scaled to 550nm. The retrieved
   exponent is used by the LUT method when no explicit override is provided.

Fallback chain: Dark Target -> Dark Pixel -> Constant (0.15)

### WVC Estimation (Continuum-Removal Band Depth)

Uses continuum-removal band depth at the 940nm and 1130nm absorption features:

1. **Widened absorption windows**: Shoulders in clean spectral regions
2. **Multi-band averaging**: 3 bands per reference point for improved SNR
3. **Air mass normalization**: Band depth / (1/cos(SZA) + 1/cos(VZA))
4. **Physically-calibrated**: scale ~28 g/cm² (940nm), ~42 g/cm² (1130nm)
5. **Weighted combination**: 940nm x0.7 normally; 1130nm x0.7 at high WVC

## References

1. Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143.

2. Hagolle, O. SMAC Python implementation. https://github.com/olivierhagolle/SMAC

3. libRadtran radiative transfer package. https://www.libradtran.org/

4. Kaufman, Y.J., et al. (1997). Operational remote sensing of tropospheric aerosol over land from EOS moderate resolution imaging spectroradiometer. *Journal of Geophysical Research*, 102(D14), 17051-17067.

5. Levy, R.C., et al. (2007). Second-generation operational algorithm: Retrieval of aerosol properties over land from inversion of Moderate Resolution Imaging Spectroradiometer spectral reflectance. *Journal of Geophysical Research*, 112, D13211.

6. Gao, B.-C. & Goetz, A.F.H. (1990). Column atmospheric water vapor and vegetation liquid water retrievals from airborne imaging spectrometer data. *Journal of Geophysical Research*, 95(D4), 3549-3564.

7. Thompson, D.R., et al. (2018). Optimal estimation for imaging spectrometer atmospheric correction. *Remote Sensing of Environment*, 216, 355-373. https://doi.org/10.1016/j.rse.2018.07.003

8. Vermote, E.F., et al. (1997). Second simulation of the satellite signal in the solar spectrum, 6S: An overview. *IEEE Transactions on Geoscience and Remote Sensing*, 35(3), 675-686.

## Authors

- **Original SMAC Algorithm**: H. Rahman and G. Dedieu (CESBIO)
- **GRASS GIS Implementation**: Yann Chemin <dr.yann.chemin@gmail.com>

## License

This is free and unencumbered software released into the public domain.

Copyright (C) 2025 GRASS Development Team

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

- [i.hyper.atcorr](https://github.com/yannforget/i.hyper.atcorr) - Companion GRASS add-on: full 6SV2.1 C port
  with BRDF surface models, custom Mie aerosol, and ISOFIT improvements.
  Use as a self-contained alternative to the libRadtran LUT backend.
- [i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html) - GRASS 6S-based atmospheric correction (multispectral)
- [i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html) - Vegetation indices
- [r3.in.xyz](https://grass.osgeo.org/grass-stable/manuals/r3.in.xyz.html) - Import 3D raster data
- [ISOFIT](https://github.com/isofit/isofit) - Imaging Spectrometer Optimal FITting
