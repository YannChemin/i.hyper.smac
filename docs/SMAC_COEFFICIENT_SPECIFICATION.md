# SMAC Coefficient File Specification

## Overview

The SMAC (Simplified Method for Atmospheric Correction) algorithm requires 49 coefficients organized in 19 lines. These coefficients are derived by fitting simple analytical formulas to radiative transfer simulations using models like 6S or libRadtran.

## References

- Rahman, H., & Dedieu, G. (1994). SMAC: a simplified method for the atmospheric correction of satellite measurements in the solar spectrum. *Int. J. Remote Sensing*, 15(1), 123-143.
- [Original SMAC implementation by Olivier Hagolle](https://github.com/olivierhagolle/SMAC)
- [libRadtran software package](https://www.libradtran.org/)

## Coefficient File Format

The coefficient file contains 19 lines with the following structure:

```
Line  0: ah2o nh2o                    # H2O absorption coefficients (2 values)
Line  1: ao3 no3                      # O3 absorption coefficients (2 values)
Line  2: ao2 no2 po2                  # O2 absorption coefficients (3 values)
Line  3: aco2 nco2 pco2               # CO2 absorption coefficients (3 values)
Line  4: ach4 nch4 pch4               # CH4 absorption coefficients (3 values)
Line  5: ano2 nno2 pno2               # NO2 absorption coefficients (3 values)
Line  6: aco nco pco                  # CO absorption coefficients (3 values)
Line  7: a0s a1s a2s a3s              # Spherical albedo coefficients (4 values)
Line  8: a0T a1T a2T a3T              # Transmission coefficients (4 values)
Line  9: taur sr                      # Rayleigh optical thickness (2 values)
Line 10: a0taup a1taup                # Aerosol optical depth coefficients (2 values)
Line 11: wo gc                        # Single scattering albedo, asymmetry (2 values)
Line 12: a0P a1P a2P                  # Phase function coefficients (3 values)
Line 13: a3P a4P                      # Phase function coefficients (2 values)
Line 14: Rest1 Rest2                  # Residual terms (2 values)
Line 15: Rest3 Rest4                  # Residual terms (2 values)
Line 16: Resr1 Resr2 Resr3            # Rayleigh residual terms (3 values)
Line 17: Resa1 Resa2                  # Aerosol residual terms (2 values)
Line 18: Resa3 Resa4                  # Aerosol residual terms (2 values)
```

**Total: 49 coefficients**

## Coefficient Meanings and Formulas

### 1. Gaseous Absorption Coefficients (Lines 0-6)

Gaseous transmission is modeled as:

```
T_gas = exp(a * (u * m)^n)
```

Where:
- `T_gas` = gaseous transmission
- `u` = column content of the gas (e.g., cm-atm for O3, g/cm² for H2O)
- `m` = air mass = 1/cos(θs) + 1/cos(θv)
- `a` = absorption coefficient (negative value)
- `n` = exponent (typically 0.5-1.0)

For pressure-dependent gases (O2, CO2, CH4, NO2, CO):
```
u_gas = P_eq^p
```
Where `P_eq = P/1013.25` is the pressure ratio and `p` is the pressure exponent.

| Gas | Coefficients | Typical Values | Absorption Bands |
|-----|--------------|----------------|------------------|
| H2O | ah2o, nh2o | ah2o < 0, nh2o ≈ 0.5 | 720, 820, 940, 1130, 1380, 1900nm |
| O3  | ao3, no3 | ao3 < 0, no3 ≈ 1.0 | 500-700nm (Chappuis band) |
| O2  | ao2, no2, po2 | Often zero for visible | 762nm (A-band) |
| CO2 | aco2, nco2, pco2 | Often zero except SWIR | 1600, 2000nm |
| CH4 | ach4, nch4, pch4 | Often zero except SWIR | 1660, 2300nm |
| NO2 | ano2, nno2, pno2 | Often zero except blue | 400-500nm |
| CO  | aco, nco, pco | Often zero except SWIR | 2300nm |

### 2. Spherical Albedo Coefficients (Line 7)

The spherical albedo of the atmosphere is modeled as:

```
s = a0s * P_eq + a3s + a1s * τ_550 + a2s * τ_550²
```

Where:
- `s` = spherical albedo
- `P_eq` = pressure ratio (P/1013.25)
- `τ_550` = aerosol optical depth at 550nm

### 3. Transmission Coefficients (Line 8)

Total scattering transmission:

```
T_down = a0T + a1T * τ_550/cos(θs) + (a2T * P_eq + a3T)/(1 + cos(θs))
T_up   = a0T + a1T * τ_550/cos(θv) + (a2T * P_eq + a3T)/(1 + cos(θv))
```

### 4. Rayleigh Optical Thickness (Line 9)

```
τ_r = taur * P_eq
```

The Rayleigh optical thickness at sea level can be computed analytically (Hansen & Travis, 1974):

```
τ_r(λ) = 0.008569 * λ^(-4) * (1 + 0.0113*λ^(-2) + 0.00013*λ^(-4))
```

Where λ is wavelength in micrometers.

| Wavelength | τ_r (sea level) |
|------------|-----------------|
| 400nm | 0.378 |
| 450nm | 0.221 |
| 550nm | 0.097 |
| 650nm | 0.049 |
| 850nm | 0.017 |
| 1000nm | 0.009 |

### 5. Aerosol Optical Depth Scaling (Line 10)

Spectral AOD from reference wavelength (550nm):

```
τ_aer(λ) = a0taup + a1taup * τ_550
```

Typically follows Ångström law:
```
τ_aer(λ) = τ_550 * (λ/550)^(-α)
```
Where α is the Ångström exponent (typically 1.0-1.5 for continental aerosols).

| Wavelength | a1taup (α=1.3) |
|------------|----------------|
| 400nm | 1.52 |
| 550nm | 1.00 |
| 850nm | 0.58 |
| 1600nm | 0.27 |

### 6. Aerosol Properties (Line 11)

- `wo` = Single scattering albedo (typically 0.85-0.99)
- `gc` = Asymmetry parameter (typically 0.6-0.8)

### 7. Phase Function Coefficients (Lines 12-13)

The aerosol phase function is approximated as a polynomial in scattering angle:

```
P(Θ) = a0P + a1P*Θ + a2P*Θ² + a3P*Θ³ + a4P*Θ⁴
```

Where Θ is the scattering angle in degrees.

For Henyey-Greenstein approximation:
```
a0P = (1 - g²) / (1 + g)^1.5
```

### 8. Residual Terms (Lines 14-18)

These coefficients account for:
- **Rest1-Rest4**: Coupling between Rayleigh and aerosol scattering
- **Resr1-Resr3**: Rayleigh scattering residuals
- **Resa1-Resa4**: Aerosol scattering residuals

These are typically small correction terms (<0.05).

## Coefficient Generation

### Method 1: Using libRadtran (Recommended)

The `lib/smac_coef_generator.py` module provides full libRadtran integration:

```python
from lib.smac_coef_generator import SMACCoefficientGenerator

generator = SMACCoefficientGenerator(verbose=True)
coef = generator.generate_coefficients(
    wavelength_nm=550.0,
    fwhm_nm=10.0,
    aerosol_type='continental'
)
coef.to_file('coef_550nm_CONTINENTAL.dat')
```

#### libRadtran Simulation Setup

The generator runs uvspec simulations with:
- Atmosphere: US Standard (afglus.dat)
- Solar spectrum: Kurucz 1.0nm resolution
- RTE solver: DISORT with 16 streams
- Output: edir, edn, eup, eglo, uu

#### Gas Absorption Fitting

For each gas, simulations vary the column content while measuring surface direct irradiance:

```python
# O3 column range: 0.2-0.45 cm-atm (200-450 DU)
# H2O column range: 0.5-5.0 g/cm²

# Fit: T = exp(a * (u*m)^n)
# where u = column content, m = air mass
```

#### Spherical Albedo Fitting

Run simulations with varying surface albedo (0.0-0.5) and AOD (0.0-0.5):

```python
# s = (E_glo(ρ₁) - E_glo(ρ₂)) / (ρ₁*E_glo(ρ₁) - ρ₂*E_glo(ρ₂))
```

### Method 2: Analytical Formulas (Fallback)

When libRadtran is unavailable, analytical formulas provide approximate coefficients:

```python
from lib.smac_coef_generator import create_analytical_coefficients

coef = create_analytical_coefficients(
    wavelength_nm=550.0,
    aerosol_type='continental'
)
```

The analytical method uses:
- Hansen & Travis (1974) formula for Rayleigh optical thickness
- Ångström law for AOD spectral dependence
- Gaussian absorption profiles for O3 (Chappuis band) and H2O bands
- Henyey-Greenstein phase function approximation

### Batch Generation

Use the `scripts/generate_hyperspectral_coefs.py` script:

```bash
# For a specific sensor
python generate_hyperspectral_coefs.py --sensor PRISMA --aerosol continental

# For custom wavelength range
python generate_hyperspectral_coefs.py --start 400 --end 2500 --step 10

# Analytical mode (no libRadtran required)
python generate_hyperspectral_coefs.py --sensor AVIRIS --analytical

# All aerosol models
python generate_hyperspectral_coefs.py --sensor HYPERION --all-aerosols
```

## Aerosol Model Types

| Model | wo | gc | Description |
|-------|----|----|-------------|
| Continental | 0.89 | 0.65 | Rural, dust mixture |
| Maritime | 0.98 | 0.72 | Sea salt dominated |
| Urban | 0.82 | 0.62 | Soot-rich |
| Desert | 0.92 | 0.72 | Mineral dust |

## Example Coefficient File

550nm, Continental aerosol (libRadtran-generated):

```
0.000000e+00 5.000000e-01
-4.561000e-02 1.003223e+00
0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00
4.863751e-02 2.000000e-01 -5.000000e-02 0.000000e+00
1.000000e+00 -1.500000e-01 0.000000e+00 -1.000000e-01
9.727502e-02 9.727502e-02
0.000000e+00 1.000000e+00
8.900000e-01 6.500000e-01
2.724746e-01 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00
```

## Pre-generated Coefficient Files

The `COEFS/` directory contains coefficient files for:

| Directory | Aerosol Type | Files | Wavelength Range |
|-----------|--------------|-------|------------------|
| CONTINENTAL/ | Continental | 43 | 400-2500nm @ 50nm |
| MARITIME/ | Maritime | 43 | 400-2500nm @ 50nm |
| URBAN/ | Urban | 43 | 400-2500nm @ 50nm |
| DESERT/ | Desert | 43 | 400-2500nm @ 50nm |

**Total: 172 coefficient files**

File naming convention: `coef_{wavelength}nm_{AEROSOL_TYPE}.dat`

## Validation

### Accuracy Targets

Generated coefficients should achieve:
- Reflectance error < 3% for typical conditions (AOD < 0.3, SZA < 60°)
- Reflectance error < 5% for moderate conditions (AOD < 0.5, SZA < 70°)

### Validation Procedure

1. Compare SMAC output with full libRadtran simulations
2. Test across range of:
   - Solar zenith angles: 0°, 30°, 45°, 60°
   - View zenith angles: 0°, 15°, 30°
   - AOD values: 0.0, 0.1, 0.2, 0.3, 0.5
   - Surface reflectances: 0.05, 0.1, 0.2, 0.3

### libRadtran vs Analytical Comparison

| Wavelength | taur (libRadtran) | taur (Analytical) | Difference |
|------------|-------------------|-------------------|------------|
| 450nm | 0.221292 | 0.220629 | 0.3% |
| 550nm | 0.097275 | 0.097142 | 0.1% |
| 650nm | 0.049323 | 0.049288 | 0.1% |
| 850nm | 0.016676 | 0.016672 | 0.02% |

The Rayleigh optical thickness from libRadtran matches analytical formulas to within 0.5%.

## Code Reference

Key files for coefficient generation:

| File | Purpose |
|------|---------|
| `lib/smac_coef_generator.py` | Main coefficient generator |
| `scripts/generate_hyperspectral_coefs.py` | Batch generation script |
| `lib/smac.py` | SMAC algorithm using coefficients |
| `tests/test_smac.py` | Unit tests for SMAC functions |
