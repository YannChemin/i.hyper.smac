# i.hyper.smac

## NAME

**i.hyper.smac** - Atmospheric correction of hyperspectral imagery using SMAC

## KEYWORDS

imagery, atmospheric correction, hyperspectral, reflectance, SMAC

## SYNOPSIS

**i.hyper.smac**

**i.hyper.smac --help**

**i.hyper.smac** input=*name* output=*name* dem=*name* solar_zenith=*float* solar_azimuth=*float* view_zenith=*float* view_azimuth=*float* [method=*string*] [sensor=*string*] [aod=*float*] [water_vapor=*float*] [ozone=*float*] [pressure=*float*] [aerosol_model=*string*] [angstrom=*float*] [visibility=*float*] [adjacency_psf=*float*] [output_uncertainty=*name*] [**-k**] [**-g**] [**-p**] [**-c**] [**-u**] [--overwrite] [--verbose] [--quiet]

## DESCRIPTION

**i.hyper.smac** performs atmospheric correction on hyperspectral imagery using the Simplified Method for Atmospheric Correction (SMAC). It converts top-of-atmosphere (TOA) radiance to surface reflectance by removing atmospheric effects.

**Current Implementation**: The module currently implements only the LUT (look-up table) method using libRadtran DISORT simulations. The simple and libradtran methods referenced in the documentation are planned for future releases.

The module corrects for:

- **Rayleigh scattering**: Molecular scattering that increases toward shorter wavelengths
- **Aerosol scattering and absorption**: Particle-related atmospheric effects
- **Gaseous absorption**: H2O, O3, O2, CO2, CH4, NO2, CO absorption bands
- **Adjacency effects**: Photon scattering from neighbouring pixels (optional)

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

### LUT Method (Recommended)

The **lut** method (default) uses a precomputed 3D look-up table from libRadtran
DISORT simulations. The LUT is indexed by (wavelength, AOD, H2O) and provides
R_atm, T_down, T_up, and s at each grid point. Trilinear interpolation delivers
per-pixel atmospheric parameters.

Key features:

- **6S-style two-albedo decomposition**: Separates downward and upward transmittance
  using albedo=0 and albedo=0.5 simulations
- **H2O as LUT dimension**: 7-point H2O grid with per-pixel water vapor retrieved
  from the 940nm absorption band (ISOFIT-style), eliminating the need for separate
  gas transmission corrections. Transmittance is interpolated in log-space for
  accurate SWIR continuum representation.
- **Path radiance and transmittance calibration**: Dark vegetation targets calibrate
  the LUT path radiance in the VIS and transmittance in the NIR/SWIR
- **Native aerosol spectral shape**: Uses Shettle/Fenn aerosol models with their
  native spectral properties, with scene-retrieved Angstrom exponent as fallback
- **LUT caching**: Generated LUTs saved as .npz files for reuse

### ISOFIT-Inspired Improvements

The LUT method includes six improvements inspired by the
[ISOFIT](https://github.com/isofit/isofit) optimal estimation framework
(Thompson et al., 2018):

1. **Superpixel atmospheric retrieval**: AOD and H2O are retrieved on SLIC
   superpixel means (high SNR) and spatially interpolated to every pixel via
   Gaussian smoothing. Reduces spatial noise in atmosphere estimates.

2. **In-loop adjacency correction**: The Vermote et al. (1997) adjacency effect
   correction is applied inside the per-band inversion loop rather than as
   post-processing. Environmental reflectance is computed via spatial averaging
   with a point spread function, and the diffuse transmittance component is
   corrected for photon scattering from neighbouring pixels.

3. **Surface reflectance prior**: A 3-component Gaussian mixture model
   (vegetation, soil, water) provides per-band prior mean and variance. Pixels
   are classified to the nearest component using VNIR bands.

4. **Uncertainty propagation**: Per-band, per-pixel reflectance uncertainty is
   computed from instrument noise (NEDL), AOD retrieval uncertainty, and
   systematic RT model errors. Output as a parallel uncertainty raster.

5. **MAP inner loop**: The surface prior and measurement uncertainty are combined
   in a maximum a posteriori (MAP) estimate. The prior has stronger influence in
   noisy bands (where measurement uncertainty is large relative to prior variance),
   providing spectral regularization.

6. **Model discrepancy matrix**: A per-band noise floor captures systematic RT
   model errors. Higher discrepancy at gas absorption band edges (720, 760, 940,
   1135nm) and in the SWIR prevents overfitting in unreliable bands.

### Automatic Parameter Estimation

When atmospheric parameters are not provided, the module estimates them
from the hyperspectral data itself:

#### AOD (Aerosol Optical Depth)

Implements the MODIS Dark Target algorithm (Kaufman 1997, Levy 2007) adapted
for hyperspectral data:

- Extracts four key bands (blue 470nm, red 660nm, NIR 860nm, SWIR 2130nm)
  and converts radiance to TOA reflectance
- Selects Dark Target pixels using NDVI and SWIR reflectance thresholds,
  filtered to the 20th-50th percentile of SWIR reflectance
- Predicts surface VIS reflectance from SWIR using empirical ratios
  (Kaufman 1997): rho_blue = 0.25 x rho_SWIR, rho_red = 0.50 x rho_SWIR
- Inverts atmospheric path radiance to AOD using single-scattering
  approximation with Henyey-Greenstein phase function
- Derives the Angstrom exponent from blue/red AOD ratio and scales to 550nm.
  The retrieved exponent is used by the LUT method when no explicit `angstrom`
  value is provided.
- Fallback chain: Dark Target -> Dark Pixel -> Constant (0.15)

#### Water Vapor Content (WVC)

Uses continuum-removal band depth at the 940nm and 1130nm water vapor
absorption features:

- Widened absorption windows with shoulder bands placed in clean spectral
  regions (940nm: 860-1010nm, 1130nm: 1060-1200nm)
- Multi-band averaging: 3 bands averaged at each reference point (left
  shoulder, absorption center, right shoulder) for improved SNR
- Air mass normalization: band depth divided by two-way air mass factor
  (1/cos(SZA) + 1/cos(VZA)) to obtain vertical column equivalent
- Physically-calibrated coefficients based on H2O absorption cross-sections
  (scale ~28 for 940nm, ~42 for 1130nm)
- Weighted combination: 940nm weighted 0.7 in normal conditions; 1130nm
  weighted 0.7 when WVC > 3.5 g/cm2 (940nm feature saturates at high WVC)

#### Ozone

Uses the Chappuis band method (500-700nm), converted from DU to cm-atm

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

**solar_azimuth**=*float*
:   Solar azimuth angle in degrees (0-360)

**view_zenith**=*float*
:   View (sensor) zenith angle in degrees (0-90)

**view_azimuth**=*float*
:   View azimuth angle in degrees (0-360)

### Optional Parameters

**method**=*string*
:   Atmospheric correction method. Options: lut (currently implemented), simple, libradtran (planned). Default: lut

**sensor**=*string*
:   Sensor type for automatic band configuration. Options: PRISMA, AVIRIS, AVIRIS_NG, HYPERION, ENMAP, OSK_GHOST, PIXXEL, ESPER, IPERLITE, KUVASPACE_23, KUVASPACE_32, WYVERN_23, WYVERN_32, HYP4U, TANAGER

**aod**=*float*
:   Aerosol optical depth at 550nm. If not specified, estimated automatically using the Dark Target algorithm

**water_vapor**=*float*
:   Water vapor column content in g/cm^2. If not specified, estimated from 940nm/1130nm absorption features with air mass correction

**ozone**=*float*
:   Ozone column content in cm-atm. Default: 0.3

**pressure**=*float*
:   Atmospheric pressure in hPa. If not specified, estimated from DEM

**aerosol_model**=*string*
:   Aerosol model type. Options: continental, maritime, urban, desert. Default: continental

**angstrom**=*float*
:   Angstrom exponent for aerosol wavelength dependence. If not specified, uses the scene-retrieved value from the Dark Target algorithm (when available), or falls back to the native Shettle/Fenn spectral shape for the selected aerosol model.

**visibility**=*float*
:   Visibility in km. If not provided, estimated from AOD

**adjacency_psf**=*float*
:   Adjacency effect PSF radius in km (0 = disabled, typical 1.0). For the LUT method, correction is applied in-loop during inversion; for other methods, applied as post-processing. Default: 0

**output_uncertainty**=*name*
:   Name for output uncertainty 3D raster map. Requires **-u** flag.

### Flags

**-k**
:   Keep temporary bands (for debugging)

**-g**
:   Generate SMAC coefficients from libRadtran at runtime instead of using pre-generated files. Requires libRadtran and scipy. Only applies to method=libradtran. Generated coefficients are cached to the COEFS directory for reuse in subsequent runs. (Planned feature)

**-p**
:   Apply spectral polishing to remove outlier bands using moving-median detection along the spectral axis. When used with the LUT method, surface prior regularization (MAP) is applied before polishing.

**-c**
:   Clear cached LUT and regenerate from libRadtran

**-u**
:   Compute per-band reflectance uncertainty and write to the **output_uncertainty** raster. Uncertainty includes instrument noise, AOD retrieval error, and model discrepancy contributions. Also enables surface prior regularization with proper measurement variance weighting.

## EXAMPLES

### Basic LUT correction with automatic parameter estimation

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0
```

### LUT correction with adjacency effect and spectral polishing

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0 adjacency_psf=1.0 -p
```

### LUT correction with uncertainty output

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0 \
    output_uncertainty=tanager_unc -u
```

### Full ISOFIT-inspired processing

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0 adjacency_psf=1.0 \
    output_uncertainty=tanager_unc -p -u
```

### Correction with specified atmospheric parameters

```bash
i.hyper.smac input=prisma_toa output=prisma_sr dem=elevation \
    solar_zenith=28 solar_azimuth=145 view_zenith=0 view_azimuth=0 \
    aod=0.12 water_vapor=1.8 ozone=0.32
```

### Maritime aerosol model for coastal scenes

```bash
i.hyper.smac input=coastal_toa output=coastal_sr dem=elevation \
    solar_zenith=40 solar_azimuth=90 view_zenith=0 view_azimuth=0 \
    aerosol_model=maritime aod=0.08 water_vapor=3.5
```

### Off-nadir viewing geometry

```bash
i.hyper.smac input=offnadir_toa output=offnadir_sr dem=elevation \
    solar_zenith=30 solar_azimuth=120 view_zenith=15 view_azimuth=270 \
    aod=0.2
```

### Force LUT regeneration

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0 -c
```

### Legacy libradtran method with runtime coefficient generation

```bash
i.hyper.smac input=tanager_toa output=tanager_sr dem=elevation \
    solar_zenith=29 solar_azimuth=180 view_zenith=0 view_azimuth=0 \
    method=libradtran aerosol_model=continental -g
```

## NOTES

### Input Data Requirements

- Input must be a 3D raster (GRASS RASTER3D) with spectral bands as the third dimension
- Values should be in TOA radiance (W/m^2/sr/um)
- Band wavelengths and FWHM should be stored in the raster metadata (as set by i.hyper.import)

### LUT Caching

The LUT method generates atmospheric look-up tables from libRadtran DISORT
simulations. Generated LUTs are cached as `.npz` files in `COEFS/LUT/` and
automatically reused when the scene geometry and atmospheric parameters match.
Use the **-c** flag to force regeneration.

Each LUT contains 9 AOD values x 7 H2O values x 2 albedos = 126 libRadtran
runs. The dense 7-point H2O grid (0.3x-2.5x scene mean) accurately captures
the nonlinear H2O continuum absorption in SWIR bands (1500-2500nm).
Transmittance parameters (T_down, T_up) are interpolated in log-space,
since ln(T) is nearly linear in H2O column, eliminating the systematic bias
from linear interpolation on the convex exp(-tau) function.
Generation takes approximately 10-20 minutes depending on wavelength range
and step size.

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

### Adjacency Effect Correction

The adjacency correction follows Vermote et al. (1997). When **adjacency_psf**
is set to a positive value (typically 1.0 km), the module computes environmental
reflectance via spatial averaging and decomposes the transmittance into direct
(Beer-Lambert) and diffuse components. The correction accounts for photon
scattering from neighbouring pixels that contributes through the diffuse
transmittance.

For the LUT method, adjacency correction is applied in-loop (per-band, during
inversion), so that subsequent spectral polishing and surface prior
regularization see the corrected values. For other methods, it is applied as
post-processing.

### Aerosol Spectral Calibration

The LUT method applies three calibration techniques to reduce systematic
reflectance bias in the NIR/SWIR:

1. **Native Shettle/Fenn spectral shape**: Aerosol models use their internally
   consistent spectral properties from libRadtran's `aerosol_haze` directive,
   without overriding the Angstrom exponent. This preserves the physically-based
   spectral variation of aerosol optical depth across the full wavelength range.

2. **Scene-retrieved Angstrom exponent**: The Dark Target AOD estimation derives
   the Angstrom exponent from the blue/red AOD ratio. This retrieved value is
   used as the LUT Angstrom parameter when the user has not provided an explicit
   override, allowing the LUT aerosol spectral shape to adapt to observed
   conditions.

3. **NIR transmittance correction**: Dark vegetation targets at 860nm (with known
   surface reflectance ~ 3.5 x rho_SWIR) are used to derive an empirical
   transmittance correction factor. The ratio of LUT-inverted to expected NIR
   surface reflectance is applied to T_down x T_up for all bands beyond 700nm,
   correcting residual aerosol transmittance bias.

### Superpixel Smoothing

The LUT method automatically performs SLIC superpixel segmentation on a NIR
reference band (~860nm) and smooths the retrieved AOD and H2O fields over
the superpixel means. This reduces spatial noise in atmospheric parameter
estimates while preserving scene-level gradients. Requires `scikit-image`.

### Surface Prior and MAP Regularization

When spectral polishing (**-p**) or uncertainty (**-u**) is enabled, the LUT
method applies a surface reflectance prior using a 3-component Gaussian mixture
model (vegetation, soil, water). Each pixel is classified to the nearest
component, and a diagonal-covariance MAP estimate blends the observation with
the prior. The prior has stronger influence in noisy bands (where measurement
uncertainty is large), providing spectral regularization without ad-hoc
thresholds.

### Uncertainty Output

When the **-u** flag and **output_uncertainty** are specified, the module writes
a parallel 3D raster containing per-band, per-pixel reflectance uncertainty
(sigma_rho). The uncertainty combines three sources in quadrature:

- **Instrument noise**: NEDL estimated from dark pixels, propagated through the inversion
- **AOD uncertainty**: +/-0.04 AOD perturbation propagated through LUT interpolation
- **Model discrepancy**: Per-band systematic RT error floor, higher at gas absorption band edges (720, 760, 940, 1135nm) and in the SWIR

### Aerosol Models

| Model | Description | Typical Conditions |
|-------|-------------|-------------------|
| Continental | Rural/background aerosol | Inland areas, moderate pollution |
| Maritime | Sea salt dominated | Coastal and oceanic areas |
| Urban | Soot and industrial particles | Cities, industrial regions |
| Desert | Mineral dust | Arid regions, dust events |

### Performance Considerations

- The **lut** method is the most accurate and is recommended for general use
- The **simple** method is fastest but least accurate
- The **-g** flag adds overhead (several seconds per band for libRadtran fitting) but generated coefficients are cached for reuse
- Automatic AOD/WVC estimation adds overhead but improves accuracy
- Superpixel smoothing requires `scikit-image`; if not installed, falls back to raw per-pixel estimates
- Uncertainty computation (**-u**) adds ~20% overhead per band
- For time-series processing, consider pre-computing atmospheric parameters

### Limitations

- Assumes horizontally homogeneous atmosphere (partially mitigated by per-pixel retrieval and adjacency correction)
- Best accuracy for AOD < 0.5 and SZA < 70 degrees
- May underperform in presence of clouds, shadows, or extreme aerosol loading
- Water bodies and snow may not be accurately corrected (surface prior defaults to water component for very dark pixels)

## REFERENCES

Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *International Journal of Remote Sensing*, 15(1), 123-143. https://doi.org/10.1080/01431169408954055

Hagolle, O. SMAC Python implementation. GitHub repository. https://github.com/olivierhagolle/SMAC

libRadtran radiative transfer package. https://www.libradtran.org/

Hansen, J.E. & Travis, L.D. (1974). Light scattering in planetary atmospheres. *Space Science Reviews*, 16, 527-610.

Kaufman, Y.J., et al. (1997). Operational remote sensing of tropospheric aerosol over land from EOS moderate resolution imaging spectroradiometer. *Journal of Geophysical Research*, 102(D14), 17051-17067.

Levy, R.C., et al. (2007). Second-generation operational algorithm: Retrieval of aerosol properties over land from inversion of MODIS spectral reflectance. *Journal of Geophysical Research*, 112, D13211.

Gao, B.-C. & Goetz, A.F.H. (1990). Column atmospheric water vapor and vegetation liquid water retrievals from airborne imaging spectrometer data. *Journal of Geophysical Research*, 95(D4), 3549-3564.

Thompson, D.R., et al. (2018). Optimal estimation for imaging spectrometer atmospheric correction. *Remote Sensing of Environment*, 216, 355-373. https://doi.org/10.1016/j.rse.2018.07.003

Vermote, E.F., et al. (1997). Second simulation of the satellite signal in the solar spectrum, 6S: An overview. *IEEE Transactions on Geoscience and Remote Sensing*, 35(3), 675-686.

## SEE ALSO

*[i.atcorr](https://grass.osgeo.org/grass-stable/manuals/i.atcorr.html)* - Atmospheric correction using 6S model

*[i.vi](https://grass.osgeo.org/grass-stable/manuals/i.vi.html)* - Vegetation indices

*[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html)* - 3D raster information

## AUTHOR

Yann Chemin, dr.yann.chemin@gmail.com

*Based on SMAC algorithm by H. Rahman and G. Dedieu (CESBIO) and Python implementation by Olivier Hagolle*

Copyright (C) 2025 GRASS Development Team
