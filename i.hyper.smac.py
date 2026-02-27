#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# MODULE:    i.hyper.smac
# AUTHOR(S): Created for hyperspectral SMAC atmospheric correction
# PURPOSE:   Apply SMAC atmospheric correction to hyperspectral imagery
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Apply SMAC atmospheric correction to hyperspectral imagery
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: atmospheric correction
# % keyword: SMAC
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option G_OPT_R3_OUTPUT
# % key: output
# % required: yes
# % description: Output atmospherically corrected 3D raster map
# % guisection: Output
# %end

# %option G_OPT_R_INPUT
# % key: dem
# % required: yes
# % description: Digital Elevation Model (DEM) in meters
# % guisection: Input
# %end

# %option
# % key: aod
# % type: double
# % required: no
# % description: Aerosol Optical Depth at 550nm (if not provided, will be estimated)
# % guisection: Atmospheric
# %end

# %option
# % key: solar_zenith
# % type: double
# % required: no
# % description: Solar zenith angle in degrees
# % guisection: Atmospheric
# %end

# %option
# % key: solar_azimuth
# % type: double
# % required: no
# % answer: 0
# % description: Solar azimuth angle in degrees
# % guisection: Atmospheric
# %end

# %option
# % key: view_zenith
# % type: double
# % required: no
# % answer: 0
# % description: View zenith angle in degrees
# % guisection: Atmospheric
# %end

# %option
# % key: view_azimuth
# % type: double
# % required: no
# % answer: 0
# % description: View azimuth angle in degrees
# % guisection: Atmospheric
# %end

# %option
# % key: water_vapor
# % type: string
# % required: no
# % description: Water vapor content (g/cm²) or method (joint,940nm,1130nm,average)
# % guisection: Atmospheric
# %end

# %option
# % key: sensor
# % type: string
# % required: no
# % options: PRISMA,AVIRIS,AVIRIS_NG,HYPERION,ENMAP,OSK_GHOST,PIXXEL,ESPER,IPERLITE,KUVASPACE_23,KUVASPACE_32,WYVERN_23,WYVERN_32,HYP4U,TANAGER
# % description: Sensor type for pre-configured parameters
# % guisection: Sensor
# %end

# %option
# % key: visibility
# % type: double
# % required: no
# % description: Visibility (km). If not provided, will be estimated from AOD.
# % guisection: Atmospheric
# %end

# %option
# % key: aerosol_model
# % type: string
# % required: no
# % options: continental,maritime,urban,desert
# % answer: continental
# % description: Aerosol model type for atmospheric correction
# % guisection: Atmospheric
# %end

# %option
# % key: angstrom
# % type: double
# % required: no
# % description: Angstrom exponent for aerosol wavelength dependence (default: auto from aerosol model)
# % guisection: Atmospheric
# %end


# %option
# % key: ozone
# % type: double
# % required: no
# % answer: 0.3
# % description: Ozone content (cm-atm)
# % guisection: Atmospheric
# %end

# %option
# % key: pressure
# % type: double
# % required: no
# % description: Atmospheric pressure (hPa)
# % guisection: Atmospheric
# %end

# %option
# % key: opencl_device
# % type: string
# % required: no
# % options: auto,gpu,cpu
# % answer: auto
# % description: OpenCL device type for GPU acceleration
# % guisection: Performance
# %end

# %option
# % key: opencl_memory
# % type: integer
# % required: no
# % answer: 1024
# % description: OpenCL memory limit in MB (0 = unlimited)
# % guisection: Performance
# %end

# %option
# % key: parallel_lut
# % type: string
# % required: no
# % options: auto,enabled,disabled
# % answer: auto
# % description: Enable parallel LUT generation using multiple uvspec processes
# % guisection: Performance
# %end

# %option
# % key: smart_lut
# % type: string
# % required: no
# % options: auto,yes,no
# % answer: auto
# % description: Use smart LUT generation with scene-specific AOD optimization
# % guisection: Performance
# %end

# %flag
# % key: k
# % description: Keep temporary bands
# % guisection: Optional
# %end


# %flag
# % key: p
# % description: Apply spectral polishing to remove outlier bands
# % guisection: Optional
# %end

# %flag
# % key: c
# % description: Clear cached LUT and regenerate from libRadtran
# % guisection: Advanced
# %end

# %option
# % key: adjacency_psf
# % type: double
# % required: no
# % answer: 0
# % description: Adjacency effect PSF radius in km (0 = disabled, typical 1.0)
# % guisection: Advanced
# %end

# %flag
# % key: u
# % description: Compute per-band reflectance uncertainty
# % guisection: Advanced
# %end

# %option G_OPT_R3_OUTPUT
# % key: output_uncertainty
# % required: no
# % description: Output uncertainty 3D raster map (requires -u flag)
# % guisection: Output
# %end

import sys
import os
import numpy as np
import multiprocessing as mp
import grass.script as gs
from grass.pygrass.raster import RasterRow
from grass.pygrass.raster.buffer import Buffer
from pathlib import Path

# Get GISBASE (GRASS installation prefix)
gisbase = os.environ.get("GISBASE")
if gisbase is None:
    # Fallback if running inside an active GRASS session
    gisbase = gs.parse_command("g.gisenv", flags="n")["GISBASE"]

lib_path = Path(gisbase) / "etc" / "i_hyper_lib"
if gs.verbosity() > 0:
    gs.message(f"Library path: {lib_path}")

if lib_path.exists():
    sys.path.insert(0, str(lib_path))
    try:
        import smac
        import radtran
        import aod
        import wvc
        import o3
        import gas_absorption  # New: enhanced gas absorption modeling
        import opencl_accelerator  # New: OpenCL GPU acceleration
        estimate_aod = aod.estimate_aod
        estimate_wvc = wvc.estimate_wvc
        get_smac_parameters = radtran.get_smac_parameters
    except ImportError as e:
        gs.fatal(f"Cannot import required modules. Make sure wvc.py and aod.py are in {lib_path}\n"
             f"Error: {e}")


def get_raster3d_info(raster3d):
    """Get information about 3D raster."""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")

def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength from band metadata by reading directly from 3D raster info."""
    try:
        # First try to get metadata directly from the 3D raster
        info = gs.read_command('r3.info', flags='h', map=raster3d)
        wavelength = None
        fwhm = None
        
        # Parse the output to find the band's metadata
        for line in info.split('\n'):
            line = line.strip()
            if line.startswith(f'Band {band_num}:'):
                try:
                    # Example line: "Band 1: 376.44000244140625 nm, FWHM: 5.389999866485596 nm"
                    parts = line.split(':', 1)[1].split(',')
                    wavelength = float(parts[0].strip().split()[0])  # Get just the number before 'nm'
                    if 'FWHM' in parts[1]:
                        fwhm = float(parts[1].split('FWHM:')[1].strip().split()[0])
                    return wavelength, fwhm, True, "nm"
                except (ValueError, IndexError) as e:
                    gs.warning(f"Error parsing band {band_num} metadata: {e}")
                    break
        
        # If we couldn't get metadata from the 3D raster, try extracting the band
        temp_band = f"tmp_meta_{os.getpid()}_{band_num}"

        try:
            # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
            
             # Extract the band
            gs.run_command('r3.to.rast', 
                          input=raster3d,
                          output=temp_band,
                          overwrite=True,
                          quiet=True)
            
            # Set the 3D region back
            gs.run_command('g.region', raster_3d=raster3d, quiet=True)
            
            # Get metadata from the extracted 2D band
            info = gs.read_command('r.info', flags='e', map=temp_band)
            for line in info.split('\n'):
                line = line.strip()
                if 'Wavelength:' in line:
                    wavelength = float(line.split('Wavelength:')[1].strip().split()[0])
                elif 'FWHM:' in line:
                    fwhm = float(line.split('FWHM:')[1].strip().split()[0])
            
            return wavelength, fwhm, True, "nm"
            
        finally:
            # Clean up temporary band
            if gs.find_file(temp_band, element='cell')['file']:
                gs.run_command('g.remove', flags='f', type='raster', name=temp_band, quiet=True)
                
    except Exception as e:
        gs.warning(f"Could not read metadata for band {band_num}: {e}")
        return None, None, True, "nm"

def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers."""
    unit = unit.lower().strip()
    
    if unit in ['nm', 'nanometer', 'nanometers']:
        return wavelength
    elif unit in ['um', 'µm', 'micrometer', 'micrometers', 'micron', 'microns']:
        return wavelength * 1000.0
    elif unit in ['m', 'meter', 'meters']:
        return wavelength * 1e9
    else:
        gs.warning(f"Unknown wavelength unit '{unit}', assuming nanometers")
        return wavelength

def get_all_band_wavelengths(raster3d):
    """Extract all band wavelengths from 3D raster."""
    info = get_raster3d_info(raster3d)
    depths = int(info['depths'])
    
    bands = []
    gs.verbose(f"Scanning {depths} bands for wavelength metadata...")
    
    for i in range(1, depths + 1):
        wavelength, fwhm, valid, unit = parse_wavelength_from_metadata(raster3d, i)
        
        if wavelength is not None:
            wavelength_nm = convert_wavelength_to_nm(wavelength, unit)
            bands.append({
                'band_num': i,
                'wavelength': wavelength_nm,
                'fwhm': fwhm if fwhm else 10,
                'valid': valid
            })
            gs.verbose(f"Band {i}: {wavelength_nm:.2f} nm")
    
    if not bands:
        gs.fatal("No wavelength metadata found. Please use data from i.hyper.import")
    
    bands.sort(key=lambda x: x['wavelength'])
    return bands

def estimate_pressure_from_dem(dem):
    """Estimate atmospheric pressure from DEM mean elevation."""
    stats = gs.parse_command('r.univar', map=dem, flags='g')
    elevation = float(stats['mean'])
    
    # Barometric formula
    pressure = 1013.25 * (1 - 0.0065 * elevation / 288.15) ** 5.255
    
    gs.message(f"Estimated pressure from DEM: {pressure:.2f} hPa (elevation: {elevation:.1f} m)")
    return pressure

def read_raster_as_array(raster_name):
    """Read a GRASS raster map into a 2D numpy array."""
    with RasterRow(raster_name) as r:
        r.open('r')
        return np.array(r)

def compute_band_transmission(coefs, sza, vza, pressure, aod550, wvc, o3):
    """Estimate total two-way transmission for a band using SMAC coefficients.

    Replicates the gaseous and scattering transmission from smac_inv
    (lines 197-208 of smac.py) without running the full inversion.

    Args:
        coefs: SMAC coefficient object (with ah2o, nh2o, ao3, ... attributes)
        sza: Solar zenith angle in degrees
        vza: View zenith angle in degrees
        pressure: Atmospheric pressure in hPa
        aod550: Aerosol optical depth at 550 nm
        wvc: Water vapor content in g/cm²
        o3: Ozone column in cm-atm

    Returns:
        Total two-way transmission (tg * T_down * T_up).  Values near zero
        indicate opaque absorption bands.
    """
    us = np.cos(np.radians(sza))
    uv = np.cos(np.radians(vza))
    Peq = pressure / 1013.25
    m = 1.0 / us + 1.0 / uv

    # Gas transmissions (same formulas as smac.py lines 197-204)
    th2o = np.exp(coefs.ah2o * ((wvc * m) ** coefs.nh2o))
    to3 = np.exp(coefs.ao3 * ((o3 * m) ** coefs.no3))

    # O2 and other gases (pressure-dependent)
    uo2 = Peq ** coefs.po2
    to2 = np.exp(coefs.ao2 * ((uo2 * m) ** coefs.no2))

    uco2 = Peq ** coefs.pco2
    tco2 = np.exp(coefs.aco2 * ((uco2 * m) ** coefs.nco2))

    uch4 = Peq ** coefs.pch4
    tch4 = np.exp(coefs.ach4 * ((uch4 * m) ** coefs.nch4))

    uno2 = Peq ** coefs.pno2
    tno2 = np.exp(coefs.ano2 * ((uno2 * m) ** coefs.nno2))

    uco = Peq ** coefs.pco
    tco = np.exp(coefs.aco * ((uco * m) ** coefs.nco))

    tg = th2o * to3 * to2 * tco2 * tch4 * tco * tno2

    # Scattering transmission (smac.py lines 207-208)
    ttetas = (coefs.a0T + coefs.a1T * aod550 / us
              + (coefs.a2T * Peq + coefs.a3T) / (1.0 + us))
    ttetav = (coefs.a0T + coefs.a1T * aod550 / uv
              + (coefs.a2T * Peq + coefs.a3T) / (1.0 + uv))

    return tg * ttetas * ttetav


def compute_gas_transmission(coefs, sza, vza, pressure, wvc, o3):
    """Compute gas-only transmission from SMAC coefficients.

    Same gas absorption formulas as compute_band_transmission() but returns
    only the gaseous transmission tg (no scattering terms).

    Args:
        coefs: SMAC coefficient object
        sza: Solar zenith angle in degrees
        vza: View zenith angle in degrees
        pressure: Atmospheric pressure in hPa (scalar or array)
        wvc: Water vapor content in g/cm² (scalar or array)
        o3: Ozone column in cm-atm (scalar or array)

    Returns:
        Gas-only transmission tg (same shape as broadcastable inputs).
    """
    us = np.cos(np.radians(sza))
    uv = np.cos(np.radians(vza))
    Peq = pressure / 1013.25
    m = 1.0 / us + 1.0 / uv

    th2o = np.exp(coefs.ah2o * ((wvc * m) ** coefs.nh2o))
    to3 = np.exp(coefs.ao3 * ((o3 * m) ** coefs.no3))

    uo2 = Peq ** coefs.po2
    to2 = np.exp(coefs.ao2 * ((uo2 * m) ** coefs.no2))
    uco2 = Peq ** coefs.pco2
    tco2 = np.exp(coefs.aco2 * ((uco2 * m) ** coefs.nco2))
    uch4 = Peq ** coefs.pch4
    tch4 = np.exp(coefs.ach4 * ((uch4 * m) ** coefs.nch4))
    uno2 = Peq ** coefs.pno2
    tno2 = np.exp(coefs.ano2 * ((uno2 * m) ** coefs.nno2))
    uco = Peq ** coefs.pco
    tco = np.exp(coefs.aco * ((uco * m) ** coefs.nco))

    return th2o * to3 * to2 * tco2 * tch4 * tco * tno2


def compute_blue_aod_taper(wavelength, aod):
    """Reduce effective AOD for SMAC to compensate for single-scattering
    overestimation of path radiance at short wavelengths and high AOD.

    Uses a saturation model: aod_eff = aod / (1 + alpha * aod)
    where alpha increases at shorter wavelengths due to stronger
    Rayleigh-aerosol coupling.  At low AOD the correction is minimal;
    at high AOD the effective AOD saturates, matching the behaviour of
    true multiple-scattering atmospheres.

    Args:
        wavelength: Band centre wavelength in nm (scalar).
        aod: Scene AOD at 550 nm (scalar or array).

    Returns:
        Effective AOD for this band (same shape as aod).
    """
    if wavelength >= 650.0:
        return aod
    # alpha: 0 at 650nm, 2.0 at ≤400nm (linear in wavelength)
    alpha = np.minimum(2.0, 2.0 * (650.0 - wavelength) / 250.0)
    return aod / (1.0 + alpha * aod)


def compute_coupling_correction(wavelength, tg, aod550, pressure,
                                 aerosol_model='continental', k=0.07):
    """Correct gas transmission for scattering-absorption coupling.

    Multiply-scattered photons travel longer paths through absorbing gas
    than assumed by the separable model r_toa = R_atm * tg + ...
    The correction amplifies absorption proportionally to the scattering
    optical depth: tg_eff = tg^(1 + k * tau_scat).

    Args:
        wavelength: Band centre wavelength in nm (scalar).
        tg: Gas-only transmission (scalar or array).
        aod550: AOD at 550 nm (scalar or array).
        pressure: Atmospheric pressure in hPa (scalar or array).
        aerosol_model: Aerosol model name (unused, reserved).
        k: Coupling strength parameter (default 0.07).

    Returns:
        Effective gas transmission tg_eff (same shape as inputs).
    """
    wl_um = wavelength / 1000.0
    tau_r = 0.008569 * wl_um**(-4) * (1 + 0.0113 * wl_um**(-2))
    tau_r = tau_r * np.asarray(pressure) / 1013.25
    tau_a = np.asarray(aod550) * (wavelength / 550.0)**(-1.3)
    tau_scat = tau_r + tau_a
    return np.asarray(tg) ** (1.0 + k * tau_scat)


def calibrate_path_radiance(input_raster, bands, atm_lut, aod,
                            solar_zenith, view_zenith, d2):
    """Derive path radiance correction from dark vegetation targets.

    Compares LUT-predicted R_atm with observed path reflectance at dark
    vegetation pixels (DDV) to correct for aerosol SSA/phase-function
    mismatch.

    Uses blue (470nm) and red (660nm) bands to fit a two-parameter
    correction model:
        R_atm_corrected(wl) = R_atm(wl) * c * (wl/550)^delta

    where c captures the overall magnitude error (SSA too high/low) and
    delta captures the spectral tilt error (Angstrom mismatch).

    Args:
        input_raster: GRASS 3D raster name (radiance).
        bands: List of band info dicts with 'band_num', 'wavelength', 'fwhm'.
        atm_lut: AtmosphericLUT instance.
        aod: Scene-mean AOD at 550nm (scalar).
        solar_zenith: Solar zenith angle (degrees).
        view_zenith: View zenith angle (degrees).
        d2: Earth-Sun distance squared.

    Returns:
        Tuple (c, delta, t_corr): correction parameters.
        c, delta: path radiance correction R_atm * c * (wl/550)^delta
        t_corr: NIR/SWIR transmittance correction factor (multiply T_down*T_up
                by t_corr for wavelengths >= 700nm). Values > 1 increase
                transmittance (reduce retrieved reflectance).
        Returns (1.0, 0.0, 1.0) if calibration fails.
    """
    try:
        import radtran
    except ImportError:
        from lib import radtran

    BLUE_TARGET = 470.0
    RED_TARGET = 660.0
    SWIR_TARGET = 2130.0
    NIR_TARGET = 860.0

    # Find nearest bands
    def find_band(target):
        best = min(bands, key=lambda b: abs(b['wavelength'] - target))
        if abs(best['wavelength'] - target) > 30:
            return None
        return best

    blue_band = find_band(BLUE_TARGET)
    red_band = find_band(RED_TARGET)
    swir_band = find_band(SWIR_TARGET)
    nir_band = find_band(NIR_TARGET)

    if not all([blue_band, red_band, swir_band, nir_band]):
        gs.warning("Cannot calibrate path radiance: missing required bands")
        return 1.0, 0.0, 1.0, 1.0

    # Extract bands to arrays
    theta_s_rad = np.radians(solar_zenith)
    cos_sza = np.cos(theta_s_rad)
    arrays = {}
    temp_maps = []

    try:
        for name, band in [('blue', blue_band), ('red', red_band),
                            ('nir', nir_band), ('swir', swir_band)]:
            band_num = band['band_num']
            wl = band['wavelength']
            fwhm = band.get('fwhm', 10.0)
            temp_name = f"tmp_calib_{name}_{os.getpid()}"
            temp_maps.append(temp_name)

            gs.run_command('g.region', t=band_num + 0.1, b=band_num,
                           quiet=True)
            gs.run_command('r3.to.rast', input=input_raster,
                           output=temp_name, overwrite=True, quiet=True)
            gs.run_command('g.rename',
                           raster=f"{temp_name}_00001,{temp_name}",
                           overwrite=True, quiet=True)

            with RasterRow(temp_name) as r:
                r.open('r')
                rad_array = np.array(r, dtype=float)

            # Convert radiance to apparent reflectance
            try:
                e0 = radtran.E0(wl, fwhm)
                if e0 is None or e0 <= 0:
                    raise ValueError("E0 invalid")
            except Exception:
                wl_m = wl * 1e-9
                hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                B = (2 * 6.626e-34 * (2.998e8)**2
                     / (wl_m**5 * (np.exp(hc_kT) - 1.0)))
                e0 = 6.794e-5 * B * 1e-6

            arrays[name] = (np.pi * rad_array * d2) / (e0 * cos_sza)

        # Reset region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        # Select dark vegetation pixels
        swir = arrays['swir']
        nir = arrays['nir']
        red = arrays['red']

        with np.errstate(divide='ignore', invalid='ignore'):
            ndvi = (nir - red) / (nir + red + 1e-6)

        valid = (np.isfinite(swir) & np.isfinite(nir)
                 & np.isfinite(red) & np.isfinite(arrays['blue']))
        dark_mask = valid & (swir > 0.01) & (swir < 0.25) & (ndvi > 0.2)

        n_dark = int(np.sum(dark_mask))
        if n_dark < 100:
            gs.warning(f"Only {n_dark} dark target pixels found, "
                       f"skipping path radiance calibration")
            return 1.0, 0.0, 1.0

        # Filter to 20th-50th percentile of SWIR
        swir_dark = swir[dark_mask]
        p20 = np.nanpercentile(swir_dark, 20)
        p50 = np.nanpercentile(swir_dark, 50)
        dark_mask &= (swir >= p20) & (swir <= p50)

        n_filtered = int(np.sum(dark_mask))
        if n_filtered < 50:
            gs.warning("Insufficient filtered dark pixels for calibration")
            return 1.0, 0.0, 1.0

        # Median observed apparent reflectance at dark pixels
        rho_swir_dark = np.nanmedian(swir[dark_mask])
        rho_blue_toa = np.nanmedian(arrays['blue'][dark_mask])
        rho_red_toa = np.nanmedian(arrays['red'][dark_mask])

        # Kaufman (1997) surface reflectance estimates
        rho_surf_blue = 0.25 * rho_swir_dark
        rho_surf_red = 0.50 * rho_swir_dark

        # LUT predictions at scene-mean AOD
        wl_blue = blue_band['wavelength']
        wl_red = red_band['wavelength']
        R_blue, Td_blue, Tu_blue, s_blue = atm_lut.interpolate(wl_blue, aod)
        R_red, Td_red, Tu_red, s_red = atm_lut.interpolate(wl_red, aod)

        # Observed path reflectance accounting for transmittance:
        # ρ_TOA = R_atm + T_down × T_up × ρ_surf / (1 - s × ρ_surf)
        # => R_atm_obs = ρ_TOA - T_down × T_up × ρ_surf / (1 - s × ρ_surf)
        surf_contrib_blue = (Td_blue * Tu_blue * rho_surf_blue
                             / (1.0 - s_blue * rho_surf_blue))
        surf_contrib_red = (Td_red * Tu_red * rho_surf_red
                            / (1.0 - s_red * rho_surf_red))

        obs_path_blue = max(rho_blue_toa - surf_contrib_blue, 0.001)
        obs_path_red = max(rho_red_toa - surf_contrib_red, 0.001)

        if R_blue < 0.001 or R_red < 0.001:
            gs.warning("LUT R_atm near zero, skipping calibration")
            return 1.0, 0.0, 1.0

        # Correction factors at blue and red
        f_blue = obs_path_blue / R_blue
        f_red = obs_path_red / R_red

        # Fit: f(wl) = c * (wl/550)^delta
        # ln(f_blue) = ln(c) + delta * ln(wl_blue/550)
        # ln(f_red)  = ln(c) + delta * ln(wl_red/550)
        ln_f_blue = np.log(f_blue)
        ln_f_red = np.log(f_red)
        ln_wl_blue = np.log(wl_blue / 550.0)
        ln_wl_red = np.log(wl_red / 550.0)

        denom = ln_wl_blue - ln_wl_red
        if abs(denom) > 0.01:
            delta = (ln_f_blue - ln_f_red) / denom
        else:
            delta = 0.0
        c = f_blue / (wl_blue / 550.0) ** delta

        # Clamp to reasonable bounds
        c = float(np.clip(c, 0.3, 1.5))
        delta = float(np.clip(delta, -2.0, 2.0))

        gs.message("Path radiance calibration from dark targets:")
        gs.message(f"  Dark pixels used: {n_filtered}")
        gs.message(f"  SWIR dark median: {rho_swir_dark:.4f}")
        gs.message(f"  Kaufman surface est: blue={rho_surf_blue:.4f}, "
                   f"red={rho_surf_red:.4f}")
        gs.message(f"  TOA apparent refl: blue={rho_blue_toa:.4f}, "
                   f"red={rho_red_toa:.4f}")
        gs.message(f"  Surface contrib (T*rho/(1-s*rho)): "
                   f"blue={surf_contrib_blue:.4f}, red={surf_contrib_red:.4f}")
        gs.message(f"  Blue ({wl_blue:.0f}nm): obs_path={obs_path_blue:.4f}, "
                   f"LUT_R_atm={R_blue:.4f}, ratio={f_blue:.3f}")
        gs.message(f"  Red ({wl_red:.0f}nm): obs_path={obs_path_red:.4f}, "
                   f"LUT_R_atm={R_red:.4f}, ratio={f_red:.3f}")
        gs.message(f"  Correction: R_atm * {c:.3f} * (wl/550)^{delta:.3f}")

        # --- P3: NIR transmittance correction from dark vegetation ---
        # Dense dark vegetation at 860nm has known surface reflectance.
        # Compare LUT-inverted reflectance with expected to derive
        # a transmittance correction factor for NIR/SWIR bands.
        t_corr = 1.0
        try:
            wl_nir = nir_band['wavelength']
            rho_nir_toa = np.nanmedian(arrays['nir'][dark_mask])

            # Expected NIR surface reflectance for dark vegetation
            # Empirical: NIR ≈ 2.0 × SWIR(2130nm) for typical vegetation
            # (Kaufman & Tanré 1992; value 3.5 was too high, causing t_corr
            # to always clamp at 0.70 and inflate r_boa by ~60%)
            rho_surf_nir_expected = np.clip(2.0 * rho_swir_dark, 0.20, 0.55)

            # LUT atmospheric parameters at 860nm
            R_nir, Td_nir, Tu_nir, s_nir = atm_lut.interpolate(wl_nir, aod)
            T_nir = float(Td_nir * Tu_nir) if np.ndim(Td_nir) == 0 else Td_nir * Tu_nir

            # Invert TOA to surface reflectance at NIR
            y_nir = (rho_nir_toa - float(R_nir)) / float(T_nir)
            rho_surf_nir_retrieved = y_nir / (1.0 + float(s_nir) * y_nir)

            if (rho_surf_nir_retrieved > 0.05 and
                    rho_surf_nir_expected > 0.05):
                # t_corr > 1 means T is too low (need to increase it)
                t_corr_raw = rho_surf_nir_retrieved / rho_surf_nir_expected
                # Clamp to reasonable range [0.7, 1.3]
                t_corr = float(np.clip(t_corr_raw, 0.7, 1.3))

                gs.message(f"  NIR transmittance correction (P3):")
                gs.message(f"    NIR TOA refl (dark veg): {rho_nir_toa:.4f}")
                gs.message(f"    LUT R_atm(NIR): {float(R_nir):.4f}, "
                           f"T(NIR): {float(T_nir):.4f}")
                gs.message(f"    Retrieved NIR surf: {rho_surf_nir_retrieved:.4f}")
                gs.message(f"    Expected NIR surf: {rho_surf_nir_expected:.4f}")
                gs.message(f"    T correction factor: {t_corr:.4f}")
            else:
                gs.warning("NIR transmittance correction: invalid reflectance values, "
                           "skipping")
        except Exception as e:
            gs.warning(f"NIR transmittance correction failed: {e}")

        return c, delta, t_corr

    except Exception as e:
        gs.warning(f"Path radiance calibration failed: {e}")
        return 1.0, 0.0, 1.0

    finally:
        for tmp in temp_maps:
            try:
                gs.run_command('g.remove', flags='f', type='raster',
                               name=tmp, quiet=True)
            except Exception:
                pass
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)


AOD_MAX = 1.5  # Upper bound for AOD validation
TRANSMISSION_THRESHOLD = 0.10  # Minimum two-way transmission
# LUT method uses stricter thresholds because the 5nm grid can't
# fully resolve sharp gas transitions at band edges
LUT_TRANSMISSION_MIN = 0.25   # Stricter gas threshold for LUT method
LUT_TSCAT_MIN = 0.15          # Minimum reliable T_scat from LUT
LUT_SPECTRAL_RATIO = 2.0     # Max T_scat variation over ±10nm before masking


def estimate_h2o_from_940nm(input_raster, bands, atm_lut, aod,
                             solar_zenith, view_zenith, d2):
    """Estimate per-pixel water vapor from the 940nm absorption band.

    ISOFIT-style retrieval: for each H2O value in the LUT grid, invert the
    observed radiance at 865nm (shoulder), 945nm (absorption), and 1040nm
    (shoulder) to surface reflectance.  The correct H2O produces a smooth
    continuum across the absorption feature (residual D ≈ 0).

    Args:
        input_raster: GRASS 3D raster name (radiance).
        bands: List of band info dicts with 'band_num', 'wavelength', 'fwhm'.
        atm_lut: AtmosphericLUT instance with h2o_grid.
        aod: Scene-mean AOD at 550nm (scalar).
        solar_zenith: Solar zenith angle (degrees).
        view_zenith: View zenith angle (degrees).
        d2: Earth-Sun distance squared.

    Returns:
        2D array of per-pixel H2O (g/cm²), or scalar h2o_ref if retrieval
        fails or LUT has no H2O dimension.
    """
    try:
        import radtran
    except ImportError:
        from lib import radtran

    if atm_lut.h2o_grid is None or atm_lut.R_atm.ndim != 3:
        gs.warning("LUT has no H2O dimension, using scene-mean H2O")
        return atm_lut.h2o_ref

    h2o_grid = atm_lut.h2o_grid

    # Target wavelengths for 940nm absorption retrieval
    SHOULDER_LO = 865.0   # Left shoulder
    ABSORPTION = 945.0     # Absorption center
    SHOULDER_HI = 1040.0   # Right shoulder

    # Find nearest bands
    def find_band(target, max_dist=30):
        best = min(bands, key=lambda b: abs(b['wavelength'] - target))
        if abs(best['wavelength'] - target) > max_dist:
            return None
        return best

    band_lo = find_band(SHOULDER_LO)
    band_abs = find_band(ABSORPTION)
    band_hi = find_band(SHOULDER_HI)

    if not all([band_lo, band_abs, band_hi]):
        gs.warning("Cannot find 865/945/1040nm bands for H2O retrieval, "
                   "using scene-mean H2O")
        return atm_lut.h2o_ref

    theta_s_rad = np.radians(solar_zenith)
    cos_sza = np.cos(theta_s_rad)
    temp_maps = []

    try:
        # Extract the 3 bands and convert to TOA reflectance
        refl_toa = {}
        for name, band in [('lo', band_lo), ('abs', band_abs),
                            ('hi', band_hi)]:
            band_num = band['band_num']
            wl = band['wavelength']
            fwhm = band.get('fwhm', 10.0)
            temp_name = f"tmp_h2o_{name}_{os.getpid()}"
            temp_maps.append(temp_name)

            gs.run_command('g.region', t=band_num + 0.1, b=band_num,
                           quiet=True)
            gs.run_command('r3.to.rast', input=input_raster,
                           output=temp_name, overwrite=True, quiet=True)
            gs.run_command('g.rename',
                           raster=f"{temp_name}_00001,{temp_name}",
                           overwrite=True, quiet=True)

            with RasterRow(temp_name) as r:
                r.open('r')
                rad_array = np.array(r, dtype=float)

            try:
                e0 = radtran.E0(wl, fwhm)
                if e0 is None or e0 <= 0:
                    raise ValueError("E0 invalid")
            except Exception:
                wl_m = wl * 1e-9
                hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                B = (2 * 6.626e-34 * (2.998e8)**2
                     / (wl_m**5 * (np.exp(hc_kT) - 1.0)))
                e0 = 6.794e-5 * B * 1e-6

            refl_toa[name] = (np.pi * rad_array * d2) / (e0 * cos_sza)

        # Reset region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        wl_lo = band_lo['wavelength']
        wl_abs = band_abs['wavelength']
        wl_hi = band_hi['wavelength']

        # Linear interpolation weight for continuum at absorption wavelength
        w_cont = (wl_abs - wl_lo) / (wl_hi - wl_lo)

        # For each H2O grid value, invert to surface reflectance and
        # compute absorption residual D = continuum(945) - rfl(945)
        ref_shape = refl_toa['lo'].shape
        D_stack = np.zeros((len(h2o_grid),) + ref_shape)

        for hi, h2o_val in enumerate(h2o_grid):
            rfl = {}
            for name, wl in [('lo', wl_lo), ('abs', wl_abs), ('hi', wl_hi)]:
                R_atm, T_down, T_up, s = atm_lut.interpolate(
                    wl, aod, h2o=h2o_val)
                # Invert: rfl = (rho_toa - R_atm) / (T_down * T_up + s * (rho_toa - R_atm))
                numerator = refl_toa[name] - R_atm
                denominator = T_down * T_up + s * numerator
                with np.errstate(divide='ignore', invalid='ignore'):
                    rfl[name] = numerator / denominator
                rfl[name] = np.nan_to_num(rfl[name], nan=0.0)

            # Continuum at absorption band = linear interp of shoulders
            continuum = rfl['lo'] + w_cont * (rfl['hi'] - rfl['lo'])
            D_stack[hi] = continuum - rfl['abs']

        # Per-pixel: find H2O that minimizes |D| by linear interpolation
        # D should go from positive (too little H2O) to negative (too much)
        # Find the zero-crossing
        h2o_map = np.full(ref_shape, atm_lut.h2o_ref)
        n_h2o = len(h2o_grid)

        if n_h2o >= 2:
            # Vectorized: for each adjacent pair, find pixels where D
            # crosses zero
            flat_shape = np.prod(ref_shape)
            D_flat = D_stack.reshape(n_h2o, flat_shape)
            h2o_flat = np.full(flat_shape, atm_lut.h2o_ref)
            assigned = np.zeros(flat_shape, dtype=bool)

            for i in range(n_h2o - 1):
                D_lo = D_flat[i]
                D_hi = D_flat[i + 1]
                # Zero crossing: signs differ or D is exactly zero
                cross = ((D_lo * D_hi) <= 0) & ~assigned
                if not np.any(cross):
                    continue
                dD = D_lo[cross] - D_hi[cross]
                with np.errstate(divide='ignore', invalid='ignore'):
                    frac = np.where(np.abs(dD) > 1e-10,
                                    D_lo[cross] / dD, 0.5)
                frac = np.clip(frac, 0.0, 1.0)
                h2o_flat[cross] = (h2o_grid[i]
                                   + frac * (h2o_grid[i + 1] - h2o_grid[i]))
                assigned[cross] = True

            # Pixels without zero crossing: use H2O with smallest |D|
            if not np.all(assigned):
                abs_D = np.abs(D_flat[:, ~assigned])
                best_idx = np.argmin(abs_D, axis=0)
                h2o_flat[~assigned] = h2o_grid[best_idx]

            h2o_map = h2o_flat.reshape(ref_shape)

        # Clamp to grid bounds
        h2o_map = np.clip(h2o_map, h2o_grid[0], h2o_grid[-1])

        gs.message(f"H2O retrieval from 940nm: "
                   f"mean={np.nanmean(h2o_map):.3f}, "
                   f"std={np.nanstd(h2o_map):.3f} g/cm²")

        return h2o_map

    except Exception as e:
        gs.warning(f"H2O retrieval failed: {e}, using scene-mean")
        return atm_lut.h2o_ref

    finally:
        for tmp in temp_maps:
            try:
                gs.run_command('g.remove', flags='f', type='raster',
                               name=tmp, quiet=True)
            except Exception:
                pass
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)


def apply_smac_correction_simple(input_raster, output_raster, bands,
                                  aod, water_vapor, ozone, pressure,
                                  solar_zenith, solar_azimuth,
                                  view_zenith, view_azimuth,
                                  keep_temp=False):
    """Optimized SMAC correction with 2D band processing."""
    
    gs.message("Applying optimized atmospheric correction with 2D band processing...")
    
    # Convert angles to radians once
    theta_s = np.radians(solar_zenith)
    theta_v = np.radians(view_zenith)
    cos_theta_s = np.cos(theta_s)
    cos_theta_v = np.cos(theta_v)
    m = 1.0 / cos_theta_s + 1.0 / cos_theta_v
    
    # Create a temporary directory
    temp_bands = []
    
    try:
        # Pre-compute all band corrections
        band_corrections = []
        for band in bands:
            wavelength = band['wavelength']

            # Atmospheric calculations
            tau_r = 0.008569 * (wavelength / 1000) ** (-4) * (1 + 0.0113 * (wavelength / 1000) ** (-2))
            tau_r *= pressure / 1013.25
            alpha = 1.3
            tau_a = aod * (wavelength / 550.0) ** (-alpha)
            tau = tau_r + tau_a
            
            # Gaseous transmission
            if 850 < wavelength < 1050:
                t_h2o = np.exp(-0.1 * water_vapor * m)
            elif 1050 < wavelength < 1250:
                t_h2o = np.exp(-0.15 * water_vapor * m)
            else:
                t_h2o = np.exp(-0.01 * water_vapor * m)
                
            t_o3 = np.exp(-0.05 * ozone * m) if 400 < wavelength < 700 else 1.0
            t_gas = t_h2o * t_o3
            t_down = np.exp(-tau / cos_theta_s)
            t_up = np.exp(-tau / cos_theta_v)
            t_total = t_down * t_up * t_gas
            rho_atm = 0.02 * tau
            
            # Enhanced correction for blue/green region (400-550nm)
            if wavelength < 550:
                # Rayleigh optical thickness approximation (varies with wavelength^-4)
                tau_r = 0.008569 * (wavelength/1000)**(-4) * (pressure/1013.25)
                
                # Aerosol optical thickness (Angstrom's law)
                alpha = 1.3  # Angstrom exponent, typical for continental aerosols
                tau_a = aod * (wavelength/550.0)**(-alpha)
                
                # Total optical thickness
                tau = tau_r + tau_a
                
                # Enhanced correction factor for blue region
                if wavelength < 450:  # Stronger correction for blue
                    enhancement = 1.0 + (450 - wavelength) * 0.01  # Adjust multiplier as needed
                else:  # Moderate correction for green
                    enhancement = 1.0 + (550 - wavelength) * 0.005
                
                # Apply enhanced correction
                t_total *= enhancement
                rho_atm *= (1.0 + (tau - 1.0) * 0.7)  # Slight adjustment to path radiance

            band_corrections.append({
                'band_num': band['band_num'],
                't_total': t_total,
                'rho_atm': rho_atm
            })
        
        # 
        # Get wavelength information from input raster
        input_info = gs.read_command('r3.info', flags='h', map=input_raster)
        wavelength_info = [line.strip() for line in input_info.split('\n') if 'Wavelength' in line]

        # Export all bands to 2D rasters
        gs.message("Exporting bands to 2D rasters...")
        for i, band in enumerate(bands):
            band_num = band['band_num']
            temp_band = f"temp_band_{band_num}_{os.getpid()}"
            temp_bands.append(temp_band)
            
            # Set region to the specific band
            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
            
            # Export the band
            gs.run_command('r3.to.rast',
                          input=input_raster,
                          output=temp_band,
                          overwrite=True,
                          quiet=True)
            
            # The output will be named temp_band_00001, rename it
            gs.run_command('g.rename',
                          raster=f"{temp_band}_00001,{temp_band}",
                          overwrite=True,
                          quiet=True)
            
            gs.percent(i, len(bands), 1)
        
        # Process each band
        gs.message("Applying atmospheric correction to bands...")
        for i, (band, corr) in enumerate(zip(bands, band_corrections)):
            band_num = band['band_num']
            temp_band = f"temp_band_{band_num}_{os.getpid()}"
            temp_band_corr = f"{temp_band}_corr"
            # Using GRASS GIS if() syntax instead of min()/max()
            expr = f"{temp_band_corr} = float(({temp_band} - {corr['rho_atm']}) / {corr['t_total']})"
            gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=True)

            # Add wavelength and FWHM to the output band
            band_wavelength = band['wavelength']
            band_comment = f"Band {band_num}: {band_wavelength} nm"
            if 'fwhm' in band:
                band_comment += f", FWHM: {band['fwhm']} nm"
            
            # Find the corresponding wavelength info for this band
            band_wl_info = next((wl for wl in wavelength_info if f"Band {band_num}:" in wl), None)
            if band_wl_info:
                # Extract the wavelength and unit from the info
                try:
                    wl_parts = band_wl_info.split(':')[1].strip().split()
                    wavelength = float(wl_parts[0])
                    unit = wl_parts[1] if len(wl_parts) > 1 else 'nm'
                    band_comment = f"Band {band_num}: {wavelength} {unit}"
                    
                    # If there's FWHM in the original info, include it
                    if 'FWHM' in band_wl_info:
                        fwhm_part = band_wl_info.split('FWHM:')[1].strip()
                        fwhm = float(fwhm_part.split()[0])
                        band_comment += f", FWHM: {fwhm} {unit}"
                except (IndexError, ValueError) as e:
                    gs.warning(f"Could not parse wavelength info for band {band_num}: {e}")
                    # Get timestamp from input raster (use r3.timestamp for 3D rasters)
            try:
                timestamp = gs.read_command('r3.timestamp', map=input_raster)
            except:
                timestamp = ""
            
            # Clean up the temporary band
            if not keep_temp:
                gs.run_command('g.remove', flags='f', type='raster', name=temp_band, quiet=True)
            
            gs.percent(i, len(bands), 1)
        
        # Combine corrected bands back into a 3D raster
        gs.message("Combining corrected bands into 3D raster...")
        # Set the 3D region back
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)
        corrected_bands = [f"temp_band_{b['band_num']}_{os.getpid()}_corr" for b in bands]
        gs.run_command('r.to.rast3', 
                      input=','.join(corrected_bands),
                      output=output_raster,
                      overwrite=True)
        # After creating the output 3D raster, transfer metadata
        try:
            # Get metadata from input raster
            input_info = gs.read_command('r3.info', flags='h', map=input_raster).strip()
            
            # Extract Band number and fwhm
            bands_info = []
            
            # Parse the input info to get all metadata
            for line in input_info.split('\n'):
                line = line.strip()
                if line.startswith('Band ') and ('nm' in line or 'um' in line):
                    bands_info.append(line)
            
            # Build the description string in the same format as Tanager import
            desc = ["Atmospheric Correction Metadata:"]
            desc.append(f"Original raster: {input_raster}")
            desc.append(f"Solar Z: {solar_zenith}°, View Z: {view_zenith}°")
            desc.append(f"AOD: {aod}, Water Vapor: {water_vapor} g/cm², Ozone: {ozone} cm-atm")
            desc.append("Measurement: Reflectance")
            desc.append("Measurement Units: -")
            
            # Add band information
            if bands_info:
                desc.append(f"Valid Bands: {len(bands_info)}")
                desc.extend(bands_info)
            
            # Set all metadata in a single r3.support call
            gs.run_command('r3.support',
                          map=output_raster,
                          title=f"SMAC corrected {input_raster}",
                          description="\n".join(desc),
                          source1="GRASS GIS i.hyper.smac module",
                          quiet=True)
            
            # Copy timestamp from input to output if available
            try:
                timestamp = gs.read_command('r3.timestamp', map=input_raster).strip()
                if timestamp:
                    gs.run_command('r3.timestamp', map=output_raster, date=timestamp)
            except:
                pass
        
        except Exception as e:
            gs.warning(f"Could not transfer all metadata to output raster: {str(e)}")
            
        gs.percent(1, 1, 1)
        gs.message("Simple atmospheric correction complete")
        
    finally:
        # Clean up temporary files
        if not keep_temp:
            gs.message("Cleaning up temporary files...")
            for temp_band in temp_bands:
                # Remove both the original and corrected bands
                for suffix in ['', '_corr']:
                    band_name = f"{temp_band}{suffix}"
                    if gs.find_file(band_name, element='cell')['file']:
                        gs.run_command('g.remove', flags='f', type='raster', 
                                     name=band_name, quiet=True)

def apply_smac_correction_libradtran(input_raster, output_raster, bands,
                                   aod, water_vapor, ozone, pressure,
                                   solar_zenith, solar_azimuth,
                                   view_zenith, view_azimuth,
                                   sensor_type, visibility=None,
                                   aerosol_model='continental', keep_temp=False,
                                   generate_coefs=False,
                                   aod_map=None, wvc_map=None,
                                   ozone_map=None, dem=None,
                                   force_regenerate=False):
    """Apply libradtran-based SMAC atmospheric correction.

    Args:
        input_raster (str): Input 3D raster name
        output_raster (str): Output 3D raster name
        bands (list): List of band information dictionaries
        aod (float): Aerosol Optical Depth at 550nm (scene-mean fallback)
        water_vapor (float): Water vapor content (g/cm²) (scene-mean fallback)
        ozone (float): Ozone content (cm-atm) (scene-mean fallback)
        pressure (float): Atmospheric pressure (hPa) (scene-mean fallback)
        solar_zenith (float): Solar zenith angle (degrees)
        solar_azimuth (float): Solar azimuth angle (degrees)
        view_zenith (float): View zenith angle (degrees)
        view_azimuth (float): View azimuth angle (degrees)
        sensor_type (str): Sensor type (e.g., 'AVIRIS', 'PRISMA')
        visibility (float, optional): Visibility in km
        aerosol_model (str): Type of aerosol model
        keep_temp (bool): Whether to keep temporary files
        generate_coefs (bool): Generate SMAC coefficients from libRadtran
        aod_map (str, optional): GRASS raster name for per-pixel AOD at 550nm
        wvc_map (str, optional): GRASS raster name for per-pixel water vapor (g/cm²)
        ozone_map (str, optional): GRASS raster name for per-pixel ozone (DU)
        dem (str, optional): GRASS raster name for DEM (m) → per-pixel pressure
    """
    
    gs.message("Applying libradtran-based SMAC atmospheric correction...")
    
    # Get acquisition date from metadata (or use current date as fallback)
    try:
        timestamp = gs.read_command('r3.timestamp', map=input_raster).strip()
        if timestamp:
            # Parse timestamp (format: "day month year")
            from datetime import datetime
            dt = datetime.strptime(timestamp.split()[0], "%d/%m/%Y")
            year, month, day = dt.year, dt.month, dt.day
        else:
            raise ValueError("No timestamp")
    except:
        # Fallback to current date
        from datetime import datetime
        dt = datetime.now()
        year, month, day = dt.year, dt.month, dt.day
        gs.warning(f"Could not get acquisition date from metadata, using current date: {year}-{month:02d}-{day:02d}")
    
    # Convert angles to radians for calculations
    theta_s_rad = np.radians(solar_zenith)
    theta_v_rad = np.radians(view_zenith)
    
    # Calculate Earth-Sun distance
    d2 = radtran.earth_sun_distance(year, month, day) ** 2
    
    # Get wavelength information from input raster
    input_info = gs.read_command('r3.info', flags='h', map=input_raster)
    wavelength_info = [line.strip() for line in input_info.split('\n') if 'Band' in line and 'nm' in line]
    
    # Load per-pixel atmospheric maps (if provided)
    if aod_map:
        aod_array = np.clip(read_raster_as_array(aod_map), 0.001, AOD_MAX)
        gs.message(f"Loaded per-pixel AOD map: {aod_map}")
    else:
        aod_array = None

    if wvc_map:
        wvc_array = np.clip(read_raster_as_array(wvc_map), 0.1, 8.0)
        gs.message(f"Loaded per-pixel WVC map: {wvc_map}")
    else:
        wvc_array = None

    if ozone_map:
        # Ozone maps are in DU, SMAC needs cm-atm
        ozone_array = np.clip(read_raster_as_array(ozone_map) * 0.001, 0.1, 0.6)
        gs.message(f"Loaded per-pixel ozone map: {ozone_map}")
    else:
        ozone_array = None

    if dem:
        elev = read_raster_as_array(dem)
        pressure_array = 1013.25 * (1 - 0.0065 * elev / 288.15) ** 5.255
        gs.message(f"Loaded DEM for per-pixel pressure: {dem}")
    else:
        pressure_array = None

    # Store temporary band names for cleanup
    temp_bands = []
    corrected_bands = []

    # Generate a blue-only LUT for bands < 650nm (uses DISORT 16-stream
    # multiple scattering instead of SMAC's Eddington two-stream which
    # overestimates path radiance at short wavelengths)
    blue_lut = None
    has_blue_bands = any(b['wavelength'] < 650 for b in bands
                         if 300 <= b['wavelength'] <= 4000)
    if has_blue_bands:
        import lut as lut_module
        phi = abs(solar_azimuth - view_azimuth)
        gs.message("Generating blue-band LUT (350-700nm) for hybrid correction...")
        blue_lut = lut_module.AtmosphericLUT.get_or_generate(
            sza=solar_zenith, vza=view_zenith, phi=phi,
            pressure=pressure, aerosol_model=aerosol_model,
            wl_min=350, wl_max=700, wl_step=2,
            force_regenerate=force_regenerate,
        )

    try:
        # Process each band
        for i, band in enumerate(bands):
            band_num = band['band_num']
            wavelength = band['wavelength']
            fwhm = band.get('fwhm', 10.0)

            # Skip bands outside the valid range for libradtran
            if wavelength < 300 or wavelength > 4000:
                gs.warning(f"Skipping band {band_num} with wavelength {wavelength} nm (outside 300-4000 nm range)")
                continue

            gs.message(f"Processing band {band_num}: {wavelength:.2f} nm, FWHM: {fwhm:.2f} nm")

            # Extract band from 3D raster
            input_band = f"tmp_input_{os.getpid()}_{band_num}"
            temp_bands.append(input_band)

            # Set the 3D region to the specific band
            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)

            # Extract the band
            gs.run_command('r3.to.rast',
                          input=input_raster,
                          output=input_band,
                          overwrite=True,
                          quiet=True)

            # The output will be named input_band_00001, rename it
            gs.run_command('g.rename',
                          raster=f"{input_band}_00001,{input_band}",
                          overwrite=True,
                          quiet=True)

            # Read raster data into numpy array
            with RasterRow(input_band) as r:
                r.open('r')
                rad_toa_band = np.array(r)

            # Get exo-atmospheric irradiance for this band
            try:
                E0_band = radtran.E0(wavelength, fwhm)
                if E0_band is None:
                    raise ValueError("E0 returned None")
            except Exception as e:
                gs.warning(f"Could not get E0 for band {band_num} at {wavelength} nm: {e}")
                gs.warning("Using approximate E0 from blackbody model")
                wl_m = wavelength * 1e-9
                hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                B = 2 * 6.626e-34 * (2.998e8)**2 / (wl_m**5 * (np.exp(hc_kT) - 1.0))
                E0_band = 6.794e-5 * B * 1e-6

            # Convert radiance to TOA reflectance
            refl_toa_band = (np.pi * rad_toa_band * d2) / (E0_band * np.cos(theta_s_rad))

            try:
                # Get SMAC coefficients (used for gas transmission in all paths)
                coefs = get_smac_parameters(
                    wavelength=wavelength,
                    fwhm=fwhm,
                    sza=solar_zenith,
                    vza=view_zenith,
                    aod_550=aod,
                    water_vapor=water_vapor,
                    ozone=ozone,
                    pressure=pressure,
                    aerosol_model=aerosol_model,
                    generate=generate_coefs,
                    verbose=gs.verbosity() > 0
                )

                # Select per-pixel or scalar atmospheric parameters
                p_pres = pressure_array if pressure_array is not None else pressure
                p_aod = aod_array if aod_array is not None else aod
                p_wvc = wvc_array if wvc_array is not None else water_vapor
                p_o3 = ozone_array if ozone_array is not None else ozone

                # Compute gas-only transmission for masking
                tg = compute_gas_transmission(
                    coefs, solar_zenith, view_zenith, p_pres, p_wvc, p_o3
                )

                # Use blue hybrid path (LUT scattering + SMAC gas) for λ < 650nm
                use_blue_hybrid = wavelength < 650 and blue_lut is not None

                if use_blue_hybrid:
                    # LUT scattering + SMAC gas path
                    if np.ndim(tg) == 0:
                        if tg < TRANSMISSION_THRESHOLD:
                            gs.verbose(
                                f"Band {band_num} ({wavelength:.1f} nm): "
                                f"gas transmission {tg:.3f} < {TRANSMISSION_THRESHOLD}, masking as NaN"
                            )
                            refl_boa_band = np.full_like(refl_toa_band, np.nan)
                        else:
                            R_atm, T_scat, s = blue_lut.interpolate(wavelength, p_aod)
                            tg_eff = compute_coupling_correction(
                                wavelength, tg, p_aod, p_pres, aerosol_model
                            )
                            numerator = refl_toa_band - R_atm * tg_eff
                            denominator = tg_eff * T_scat + numerator * s
                            with np.errstate(divide='ignore', invalid='ignore'):
                                refl_boa_band = numerator / denominator
                            refl_boa_band = np.clip(
                                np.nan_to_num(refl_boa_band, nan=0.0), -0.01, 1.5
                            )
                    else:
                        opaque_mask = tg < TRANSMISSION_THRESHOLD
                        if np.all(opaque_mask):
                            refl_boa_band = np.full_like(refl_toa_band, np.nan)
                        else:
                            R_atm, T_scat, s = blue_lut.interpolate(wavelength, p_aod)
                            tg_eff = compute_coupling_correction(
                                wavelength, tg, p_aod, p_pres, aerosol_model
                            )
                            numerator = refl_toa_band - R_atm * tg_eff
                            denominator = tg_eff * T_scat + numerator * s
                            with np.errstate(divide='ignore', invalid='ignore'):
                                refl_boa_band = numerator / denominator
                            refl_boa_band = np.clip(
                                np.nan_to_num(refl_boa_band, nan=0.0), -0.01, 1.5
                            )
                            refl_boa_band[opaque_mask] = np.nan
                else:
                    # Original SMAC path for λ >= 650nm
                    band_T = compute_band_transmission(
                        coefs, solar_zenith, view_zenith,
                        p_pres, p_aod, p_wvc, p_o3
                    )

                    if np.ndim(band_T) == 0:
                        if band_T < TRANSMISSION_THRESHOLD:
                            gs.verbose(
                                f"Band {band_num} ({wavelength:.1f} nm): "
                                f"transmission {band_T:.3f} < {TRANSMISSION_THRESHOLD}, masking as NaN"
                            )
                            refl_boa_band = np.full_like(refl_toa_band, np.nan)
                        else:
                            refl_boa_band = smac.smac_inv(
                                r_toa=refl_toa_band,
                                tetas=solar_zenith,
                                phis=solar_azimuth,
                                tetav=view_zenith,
                                phiv=view_azimuth,
                                pressure=p_pres,
                                taup550=p_aod,
                                uo3=p_o3,
                                uh2o=p_wvc,
                                coef=coefs
                            )
                            refl_boa_band = np.clip(refl_boa_band, -0.01, 1.5)
                    else:
                        opaque_mask = band_T < TRANSMISSION_THRESHOLD
                        if np.all(opaque_mask):
                            refl_boa_band = np.full_like(refl_toa_band, np.nan)
                        else:
                            refl_boa_band = smac.smac_inv(
                                r_toa=refl_toa_band,
                                tetas=solar_zenith,
                                phis=solar_azimuth,
                                tetav=view_zenith,
                                phiv=view_azimuth,
                                pressure=p_pres,
                                taup550=p_aod,
                                uo3=p_o3,
                                uh2o=p_wvc,
                                coef=coefs
                            )
                            refl_boa_band = np.clip(refl_boa_band, -0.01, 1.5)
                            refl_boa_band[opaque_mask] = np.nan

                # Write corrected band back to a raster
                output_band = f"tmp_corr_{os.getpid()}_{band_num}"
                corrected_bands.append(output_band)
                
                # Create output raster
                ncols = refl_boa_band.shape[1] if refl_boa_band.ndim > 1 else refl_boa_band.shape[0]
                with RasterRow(output_band, mode='w', mtype='DCELL', overwrite=True) as r:
                    for row_idx, row_data in enumerate(refl_boa_band):
                        buf = Buffer((ncols,), mtype='DCELL')
                        buf[:] = row_data
                        r.put_row(buf)
                
                # Add wavelength metadata to corrected band
                band_comment = f"Band {band_num}: {wavelength:.2f} nm"
                if fwhm:
                    band_comment += f", FWHM: {fwhm:.2f} nm"
                
                gs.run_command('r.support',
                              map=output_band,
                              title=band_comment,
                              units="reflectance",
                              quiet=True)
                
                gs.percent(i, len(bands), 1)
                
            except Exception as e:
                gs.fatal(f"Error processing band {band_num}: {str(e)}")
        
        if not corrected_bands:
            gs.fatal("No bands were successfully processed")
        
        # Restore 3D region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)
        
        # Combine corrected bands back into a 3D raster
        gs.message("Combining corrected bands into 3D raster...")
        gs.run_command('r.to.rast3',
                      input=','.join(corrected_bands),
                      output=output_raster,
                      overwrite=True)
        
        # Transfer metadata to output raster
        try:
            # Build metadata description
            desc = ["Atmospheric Correction Metadata (libRadtran):"]
            desc.append(f"Original raster: {input_raster}")
            desc.append("Method: libRadtran SMAC")
            desc.append(f"Solar Z: {solar_zenith}°, View Z: {view_zenith}°")
            desc.append(f"AOD: {aod}, Water Vapor: {water_vapor} g/cm²")
            desc.append(f"Ozone: {ozone} cm-atm, Pressure: {pressure} hPa")
            desc.append(f"Aerosol model: {aerosol_model}")
            desc.append("Measurement: Reflectance (Bottom of Atmosphere)")
            desc.append("Measurement Units: -")
            
            # Add band information
            if wavelength_info:
                desc.append(f"Valid Bands: {len(corrected_bands)}")
                for i, band in enumerate(bands):
                    if i < len(corrected_bands):
                        wl_line = next((w for w in wavelength_info if f"Band {band['band_num']}:" in w), None)
                        if wl_line:
                            desc.append(wl_line)
            
            # Set metadata
            gs.run_command('r3.support',
                          map=output_raster,
                          title=f"SMAC corrected (libRadtran) {input_raster}",
                          description="\n".join(desc),
                          source1="GRASS GIS i.hyper.smac module (libRadtran)",
                          quiet=True)
            
            # Copy timestamp
            try:
                timestamp = gs.read_command('r3.timestamp', map=input_raster).strip()
                if timestamp:
                    gs.run_command('r3.timestamp', map=output_raster, date=timestamp)
            except:
                pass
                
        except Exception as e:
            gs.warning(f"Could not transfer all metadata to output raster: {str(e)}")
        
        gs.percent(1, 1, 1)
        gs.message(f"Libradtran-based atmospheric correction complete: {output_raster}")
        
    except Exception as e:
        gs.fatal(f"Error in libradtran processing: {str(e)}")
        
    finally:
        # Restore 3D region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)
        
        # Clean up temporary files
        if not keep_temp:
            gs.message("Cleaning up temporary files...")
            
            # Remove input temporary bands
            for temp_band in temp_bands:
                if gs.find_file(temp_band, element='cell')['file']:
                    gs.run_command('g.remove', flags='f', type='raster',
                                  name=temp_band, quiet=True)
            
            # Remove corrected temporary bands
            for corr_band in corrected_bands:
                if gs.find_file(corr_band, element='cell')['file']:
                    gs.run_command('g.remove', flags='f', type='raster',
                                  name=corr_band, quiet=True)
    
def apply_lut_correction(input_raster, output_raster, bands,
                         aod, water_vapor, ozone, pressure,
                         solar_zenith, solar_azimuth,
                         view_zenith, view_azimuth,
                         aerosol_model='continental', keep_temp=False,
                         aod_map=None, wvc_map=None,
                         ozone_map=None, dem=None,
                         polish=False, angstrom_alpha=None,
                         force_regenerate=False,
                         adjacency_psf_km=0.0, pixel_size=None,
                         compute_uncertainty=False,
                         output_uncertainty=None,
                         opencl_device='auto', opencl_memory=1024,
                         parallel_lut='auto', smart_lut='auto'):
    """Apply atmospheric correction using libRadtran LUT.

    Uses full multiple-scattering from libRadtran DISORT via a precomputed
    look-up table for R_atm, T_down, T_up, and s.  H2O is a LUT dimension:
    per-pixel water vapor is retrieved from the 940nm absorption band
    (ISOFIT-style) and the LUT is interpolated in (wavelength, AOD, H2O)
    space — no separate gas transmission correction needed.

    Includes ISOFIT-inspired improvements:
    1. Superpixel smoothing of AOD/H2O fields
    2. In-loop adjacency effect correction
    3. Surface reflectance prior (Gaussian mixture)
    4. Per-band reflectance uncertainty propagation
    5. MAP inner loop (prior-weighted inversion)
    6. Model discrepancy noise floor

    Args:
        input_raster (str): Input 3D raster name
        output_raster (str): Output 3D raster name
        bands (list): List of band information dictionaries
        aod (float): Aerosol Optical Depth at 550nm (scene-mean fallback)
        water_vapor (float): Water vapor content (g/cm²) (scene-mean fallback)
        ozone (float): Ozone content (cm-atm) (scene-mean fallback)
        pressure (float): Atmospheric pressure (hPa) (scene-mean fallback)
        solar_zenith (float): Solar zenith angle (degrees)
        solar_azimuth (float): Solar azimuth angle (degrees)
        view_zenith (float): View zenith angle (degrees)
        view_azimuth (float): View azimuth angle (degrees)
        aerosol_model (str): Type of aerosol model
        keep_temp (bool): Whether to keep temporary files
        aod_map (str, optional): GRASS raster name for per-pixel AOD at 550nm
        wvc_map (str, optional): GRASS raster name for per-pixel water vapor (g/cm²)
        ozone_map (str, optional): GRASS raster name for per-pixel ozone (DU)
        dem (str, optional): GRASS raster name for DEM (m) -> per-pixel pressure
        adjacency_psf_km (float): PSF radius for adjacency correction (0 = disabled)
        pixel_size (float, optional): Pixel size in meters (auto-detect from GRASS)
        compute_uncertainty (bool): Whether to compute per-band uncertainty
        output_uncertainty (str, optional): Name for uncertainty 3D raster
        opencl_device (str): OpenCL device type ('auto', 'gpu', 'cpu')
        opencl_memory (int): OpenCL memory limit in MB (0 = unlimited)
        parallel_lut (str): Parallel LUT generation ('auto', 'enabled', 'disabled')
        smart_lut (str): Smart LUT generation with scene-specific AOD ('auto', 'yes', 'no')
    """
    import lut as lut_module
    import opencl_accelerator

    gs.message("Applying libRadtran LUT atmospheric correction...")
    
    # Initialize OpenCL accelerator
    accelerator = opencl_accelerator.create_opencl_accelerator(
        device_type=opencl_device, 
        verbose=gs.verbosity() > 1
    )
    
    if accelerator.is_available():
        device_info = accelerator.get_device_info()
        gs.message(f"Using OpenCL acceleration: {device_info['name']} ({device_info['type']})")
        gs.message(f"Device: {device_info['max_compute_units']} compute units, "
                  f"{device_info['global_memory'] // (1024*1024)} MB memory")
    else:
        gs.message("OpenCL not available, using CPU processing")
        accelerator = None

    # Get acquisition date from metadata
    try:
        timestamp = gs.read_command('r3.timestamp', map=input_raster).strip()
        if timestamp:
            from datetime import datetime
            dt = datetime.strptime(timestamp.split()[0], "%d/%m/%Y")
            year, month, day = dt.year, dt.month, dt.day
        else:
            raise ValueError("No timestamp")
    except Exception:
        from datetime import datetime
        dt = datetime.now()
        year, month, day = dt.year, dt.month, dt.day
        gs.warning(f"Could not get acquisition date, using current date: "
                   f"{year}-{month:02d}-{day:02d}")

    # Calculate Earth-Sun distance
    d2 = radtran.earth_sun_distance(year, month, day) ** 2

    # Relative azimuth
    phi = abs(solar_azimuth - view_azimuth)

    # Determine wavelength range from bands
    wl_min = max(350, int(bands[0]['wavelength']) - 20)
    wl_max = min(2500, int(bands[-1]['wavelength']) + 20)

    # Configure parallel LUT generation with 75% resource usage
    if parallel_lut == 'disabled':
        # Disable parallel processing
        parallel = False
        max_workers = None
    else:
        # Enable parallel processing with 75% resource usage
        parallel = True
        if parallel_lut == 'auto' and mp.cpu_count() <= 2:
            # For systems with 2 or fewer cores, use sequential
            parallel = False
            max_workers = None
        else:
            # Use new 75% resource configuration
            try:
                import parallel_lut as parallel_module
                generator = parallel_module.create_optimal_parallel_generator(
                    resource_usage=0.75,  # 75% of resources
                    device_type='auto'    # GPU first, CPU fallback
                )
                max_workers = generator.max_workers
            except ImportError:
                # Fallback to old method
                max_workers = mp.cpu_count() - 1
    
    # Generate or load atmospheric LUT
    smart_lut_option = smart_lut
    
    if smart_lut_option in ['yes', 'auto'] and aod is not None:
        # Try to use smart LUT system
        try:
            import smart_lut
            gs.message("Using smart LUT generation with scene-specific AOD optimization")
            
            # Generate smart LUT optimized for scene AOD
            atm_lut = smart_lut.get_smart_lut_or_generate(
                scene_aod=aod,
                precision_threshold=0.05,  # 5% precision requirement
                sza=solar_zenith, vza=view_zenith, phi=phi,
                pressure=pressure, aerosol_model=aerosol_model,
                wl_min=wl_min, wl_max=wl_max, wl_step=2,
                h2o=water_vapor, o3=ozone,
                angstrom_alpha=angstrom_alpha,
                force_regenerate=force_regenerate,
                parallel=parallel,
                max_workers=max_workers
            )
            gs.message(f"Smart LUT optimized for scene AOD {aod:.3f}")
            
        except ImportError:
            if smart_lut_option == 'yes':
                gs.warning("Smart LUT module not available, falling back to standard LUT generation")
            # Fall back to standard LUT generation
            atm_lut = lut_module.AtmosphericLUT.get_or_generate(
                sza=solar_zenith, vza=view_zenith, phi=phi,
                pressure=pressure, aerosol_model=aerosol_model,
                wl_min=wl_min, wl_max=wl_max, wl_step=2,
                h2o=water_vapor, o3=ozone,
                angstrom_alpha=angstrom_alpha,
                force_regenerate=force_regenerate,
                parallel=parallel,
                max_workers=max_workers
            )
    else:
        # Use standard LUT generation
        if smart_lut_option == 'yes':
            gs.warning("Smart LUT disabled: AOD not available or explicitly disabled")
        
        atm_lut = lut_module.AtmosphericLUT.get_or_generate(
            sza=solar_zenith, vza=view_zenith, phi=phi,
            pressure=pressure, aerosol_model=aerosol_model,
            wl_min=wl_min, wl_max=wl_max, wl_step=2,
            h2o=water_vapor, o3=ozone,
            angstrom_alpha=angstrom_alpha,
            force_regenerate=force_regenerate,
            parallel=parallel,
            max_workers=max_workers
        )

    # Get wavelength information from input raster
    input_info = gs.read_command('r3.info', flags='h', map=input_raster)
    wavelength_info = [line.strip() for line in input_info.split('\n')
                       if 'Band' in line and 'nm' in line]

    # Load per-pixel atmospheric maps (if provided)
    # NaN pixels (e.g. outside Dark Target coverage) are filled with scene-mean
    if aod_map:
        aod_array = np.clip(read_raster_as_array(aod_map), 0.001, AOD_MAX)
        nan_mask = np.isnan(aod_array)
        if np.any(nan_mask):
            n_nan = int(np.sum(nan_mask))
            gs.message(f"Filling {n_nan}/{aod_array.size} NaN AOD pixels "
                       f"with scene-mean {aod:.4f}")
            aod_array[nan_mask] = aod
        gs.message(f"Loaded per-pixel AOD map: {aod_map}")
    else:
        aod_array = None

    if wvc_map:
        wvc_array = np.clip(read_raster_as_array(wvc_map), 0.1, 8.0)
        nan_mask = np.isnan(wvc_array)
        if np.any(nan_mask):
            n_nan = int(np.sum(nan_mask))
            gs.message(f"Filling {n_nan}/{wvc_array.size} NaN WVC pixels "
                       f"with scene-mean {water_vapor:.4f}")
            wvc_array[nan_mask] = water_vapor
        gs.message(f"Loaded per-pixel WVC map: {wvc_map}")
    else:
        wvc_array = None

    if ozone_map:
        ozone_array = np.clip(read_raster_as_array(ozone_map) * 0.001, 0.1, 0.6)
        nan_mask = np.isnan(ozone_array)
        if np.any(nan_mask):
            n_nan = int(np.sum(nan_mask))
            gs.message(f"Filling {n_nan}/{ozone_array.size} NaN ozone pixels "
                       f"with scene-mean {ozone:.4f}")
            ozone_array[nan_mask] = ozone
        gs.message(f"Loaded per-pixel ozone map: {ozone_map}")
    else:
        ozone_array = None

    if dem:
        elev = read_raster_as_array(dem)
        pressure_array = 1013.25 * (1 - 0.0065 * elev / 288.15) ** 5.255
        nan_mask = np.isnan(pressure_array)
        if np.any(nan_mask):
            n_nan = int(np.sum(nan_mask))
            gs.message(f"Filling {n_nan}/{pressure_array.size} NaN pressure "
                       f"pixels with {pressure:.1f} hPa")
            pressure_array[nan_mask] = pressure
        gs.message(f"Loaded DEM for per-pixel pressure: {dem}")
    else:
        pressure_array = None

    # Calibrate path radiance from dark vegetation targets
    path_corr_c, path_corr_delta, path_corr_t = calibrate_path_radiance(
        input_raster, bands, atm_lut, aod,
        solar_zenith, view_zenith, d2,
    )

    # Retrieve per-pixel H2O from 940nm absorption band
    h2o_map = estimate_h2o_from_940nm(
        input_raster, bands, atm_lut, aod,
        solar_zenith, view_zenith, d2,
    )

    # --- Improvement #1: Superpixel smoothing of AOD/H2O fields ---
    try:
        import superpixel as superpixel_module
        has_superpixel = True
    except ImportError:
        has_superpixel = False

    if has_superpixel and (
        (isinstance(h2o_map, np.ndarray) and h2o_map.ndim == 2) or
        (aod_array is not None)
    ):
        # Extract a NIR reference band (~860nm) for SLIC segmentation
        nir_band = None
        nir_band_info = None
        for b in bands:
            if 840 <= b['wavelength'] <= 880:
                nir_band_info = b
                break
        if nir_band_info is None:
            # Fallback: pick closest to 860
            nir_band_info = min(bands, key=lambda b: abs(b['wavelength'] - 860))

        nir_name = f"tmp_nir_seg_{os.getpid()}"
        try:
            gs.run_command('g.region', t=nir_band_info['band_num'] + 0.1,
                          b=nir_band_info['band_num'], quiet=True)
            gs.run_command('r3.to.rast', input=input_raster,
                          output=nir_name, overwrite=True, quiet=True)
            gs.run_command('g.rename',
                          raster=f"{nir_name}_00001,{nir_name}",
                          overwrite=True, quiet=True)
            with RasterRow(nir_name) as r:
                r.open('r')
                nir_band = np.array(r)
            gs.run_command('g.remove', flags='f', type='raster',
                          name=nir_name, quiet=True)
        except Exception as e:
            gs.warning(f"Could not extract NIR band for superpixel: {e}")

        if nir_band is not None:
            gs.message("Superpixel smoothing of atmospheric fields...")
            try:
                sp_labels = superpixel_module.segment_image(nir_band, n_segments=200)

                # Smooth H2O map
                if isinstance(h2o_map, np.ndarray) and h2o_map.ndim == 2:
                    sp_h2o = superpixel_module.superpixel_means(h2o_map, sp_labels)
                    h2o_map = superpixel_module.interpolate_superpixel_field(
                        sp_labels, sp_h2o, smoothing_sigma=2.0)
                    gs.message("  H2O field smoothed via superpixels")

                # Smooth AOD map
                if aod_array is not None:
                    sp_aod = superpixel_module.superpixel_means(aod_array, sp_labels)
                    aod_array = superpixel_module.interpolate_superpixel_field(
                        sp_labels, sp_aod, smoothing_sigma=2.0)
                    aod_array = np.clip(aod_array, 0.001, AOD_MAX)
                    gs.message("  AOD field smoothed via superpixels")
            except Exception as e:
                gs.warning(f"Superpixel smoothing failed: {e}")

    # --- Auto-detect pixel size for adjacency correction ---
    if adjacency_psf_km > 0 and pixel_size is None:
        try:
            region = gs.region()
            pixel_size = (region['ewres'] + region['nsres']) / 2.0
        except Exception:
            pixel_size = 30.0  # fallback to 30m
            gs.warning(f"Could not detect pixel size, using {pixel_size}m")

    # Store temporary band names for cleanup
    temp_bands = []
    corrected_bands = []
    # In-memory storage for spectral polishing
    band_data_list = []  # list of (band_info, refl_boa_band) tuples
    # Parallel storage for uncertainty
    uncertainty_list = []  # list of (band_info, sigma_rfl) tuples

    theta_s_rad = np.radians(solar_zenith)

    # GPU-accelerated processing function
    def process_bands_gpu(bands_subset, accelerator, atm_lut, aod_array, h2o_map, 
                          path_corr_c, path_corr_delta, path_corr_t):
        """Process multiple bands using GPU acceleration."""
        if not accelerator or not accelerator.is_available():
            return None
        
        try:
            # Collect all band data for GPU processing
            band_data = []
            for band in bands_subset:
                band_num = band['band_num']
                wavelength = band['wavelength']
                fwhm = band.get('fwhm', 10.0)
                
                # Extract band from 3D raster
                input_band = f"tmp_input_{os.getpid()}_{band_num}"
                gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
                gs.run_command('r3.to.rast', input=input_raster,
                              output=input_band, overwrite=True, quiet=True)
                gs.run_command('g.rename',
                              raster=f"{input_band}_00001,{input_band}",
                              overwrite=True, quiet=True)
                
                # Read raster data
                with RasterRow(input_band) as r:
                    r.open('r')
                    rad_toa_band = np.array(r)
                
                # Get exo-atmospheric irradiance
                try:
                    E0_band = radtran.E0(wavelength, fwhm)
                    if E0_band is None:
                        raise ValueError("E0 returned None")
                except Exception:
                    wl_m = wavelength * 1e-9
                    hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                    B = 2 * 6.626e-34 * (2.998e8)**2 / (wl_m**5 * (np.exp(hc_kT) - 1.0))
                    E0_band = 6.794e-5 * B * 1e-6
                
                # Convert radiance to TOA reflectance
                refl_toa = (np.pi * rad_toa_band * d2) / (E0_band * np.cos(theta_s_rad))
                nodata = np.isnan(rad_toa_band) | (rad_toa_band <= 0)
                refl_toa[nodata] = np.nan
                
                band_data.append({
                    'band_num': band_num,
                    'wavelength': wavelength,
                    'fwhm': fwhm,
                    'refl_toa': refl_toa,
                    'input_band': input_band
                })
            
            # Stack all bands for GPU processing
            all_refl_toa = np.stack([b['refl_toa'] for b in band_data])
            wavelengths = np.array([b['wavelength'] for b in band_data])
            
            # Prepare LUT data for GPU (simplified)
            lut_data = {
                'r_atm': np.zeros((1, len(bands_subset))),
                't_down': np.zeros((1, len(bands_subset))),
                't_up': np.zeros((1, len(bands_subset))),
                's': np.zeros((1, len(bands_subset)))
            }
            
            # Prepare atmospheric maps
            atmospheric_maps = {
                'aod': aod_array if aod_array is not None else np.array([aod]),
                'wvc': h2o_map if h2o_map is not None else np.array([water_vapor]),
                'o3': np.array([ozone]),
                'pressure': np.array([pressure])
            }
            
            # Prepare geometry
            geometry = {
                'solar_zenith': solar_zenith,
                'view_zenith': view_zenith,
                'solar_azimuth': solar_azimuth,
                'view_azimuth': view_azimuth
            }
            
            # Prepare bands info
            bands_info = [{'wavelength': b['wavelength']} for b in band_data]
            
            # Apply GPU atmospheric correction
            gpu_result = accelerator.atmospheric_correction_gpu(
                all_refl_toa, lut_data, atmospheric_maps, geometry, bands_info, 
                use_maps=True
            )
            
            if gpu_result is not None:
                # Unstack results and return band data
                results = []
                for i, band in enumerate(band_data):
                    results.append({
                        **band,
                        'refl_boa': gpu_result[i],
                        'gpu_processed': True
                    })
                return results
            else:
                return None
                
        except Exception as e:
            if gs.verbosity() > 1:
                gs.warning(f"GPU processing failed: {e}")
            return None

    try:
        for i, band in enumerate(bands):
            band_num = band['band_num']
            wavelength = band['wavelength']
            fwhm = band.get('fwhm', 10.0)

            if wavelength < 300 or wavelength > 4000:
                gs.warning(f"Skipping band {band_num} at {wavelength} nm "
                           f"(outside 300-4000 nm range)")
                continue

            gs.message(f"Processing band {band_num}: {wavelength:.2f} nm, "
                       f"FWHM: {fwhm:.2f} nm")

            # Try GPU processing for batches of bands
            if accelerator and accelerator.is_available() and i % 10 == 0:
                # Process this band and next few bands with GPU
                batch_size = min(10, len(bands) - i)
                batch_bands = bands[i:i+batch_size]
                
                gpu_results = process_bands_gpu(
                    batch_bands, accelerator, atm_lut, aod_array, h2o_map,
                    path_corr_c, path_corr_delta, path_corr_t
                )
                
                if gpu_results:
                    # GPU processing succeeded, use results and skip CPU processing
                    processing_method = "GPU" if accelerator and accelerator.is_available() else "CPU"
                    gs.message(f"{processing_method} processed {len(gpu_results)} bands")
                    
                    # Skip ahead for the processed bands
                    i += len(gpu_results) - 1
                    continue
                else:
                    # GPU processing failed, fall back to CPU
                    if accelerator and accelerator.is_available():
                        gs.verbose(f"GPU processing failed, using CPU for {len(batch_bands)} bands")
                    else:
                        gs.verbose(f"OpenCL not available, using CPU for {len(batch_bands)} bands")

            # Extract band from 3D raster
            input_band = f"tmp_input_{os.getpid()}_{band_num}"
            temp_bands.append(input_band)

            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
            gs.run_command('r3.to.rast', input=input_raster,
                          output=input_band, overwrite=True, quiet=True)
            gs.run_command('g.rename',
                          raster=f"{input_band}_00001,{input_band}",
                          overwrite=True, quiet=True)

            # Read raster data
            with RasterRow(input_band) as r:
                r.open('r')
                rad_toa_band = np.array(r)

            # Get exo-atmospheric irradiance
            try:
                E0_band = radtran.E0(wavelength, fwhm)
                if E0_band is None:
                    raise ValueError("E0 returned None")
            except Exception as e:
                gs.warning(f"Could not get E0 for band {band_num}: {e}")
                gs.warning("Using approximate E0 from blackbody model")
                wl_m = wavelength * 1e-9
                hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                B = 2 * 6.626e-34 * (2.998e8)**2 / (wl_m**5 * (np.exp(hc_kT) - 1.0))
                E0_band = 6.794e-5 * B * 1e-6

            # Convert radiance to TOA reflectance
            refl_toa = (np.pi * rad_toa_band * d2) / (E0_band * np.cos(theta_s_rad))

            # Mark nodata pixels: GRASS null (NaN) or zero radiance
            # Set to NaN so they propagate through inversion naturally
            nodata = np.isnan(rad_toa_band) | (rad_toa_band <= 0)
            refl_toa[nodata] = np.nan

            # Select per-pixel or scalar AOD for LUT interpolation
            p_aod = aod_array if aod_array is not None else aod

            try:
                # Get SMAC coefficients for gas transmission
                coefs = get_smac_parameters(
                    wavelength=wavelength, fwhm=fwhm,
                    sza=solar_zenith, vza=view_zenith,
                    aod_550=aod, water_vapor=water_vapor,
                    ozone=ozone, pressure=pressure,
                    aerosol_model=aerosol_model,
                    verbose=gs.verbosity() > 0
                )

                # Compute gas transmission at scene-mean for opacity check
                tg_ref = compute_gas_transmission(
                    coefs, solar_zenith, view_zenith, pressure,
                    water_vapor, ozone
                )

                # Mask opaque/unreliable gas bands (LUT-specific threshold)
                if tg_ref < LUT_TRANSMISSION_MIN:
                    gs.verbose(
                        f"Band {band_num} ({wavelength:.1f} nm): "
                        f"gas transmission {tg_ref:.3f} < {LUT_TRANSMISSION_MIN}, "
                        f"masking as NaN"
                    )
                    refl_boa_band = np.full_like(refl_toa, np.nan)
                    sigma_rfl_band = np.full_like(refl_toa, np.nan)
                else:
                    # Apply AOD taper in the blue to compensate for aerosol
                    # model overestimating path radiance at short wavelengths
                    p_aod_eff = compute_blue_aod_taper(wavelength, p_aod)

                    # LUT R_atm, T_down, T_up, s with per-pixel H2O
                    p_h2o = h2o_map if h2o_map is not None else water_vapor
                    R_atm, T_down, T_up, s = atm_lut.interpolate(
                        wavelength, p_aod_eff, h2o=p_h2o)

                    # Apply path radiance calibration from dark targets
                    # (VIS only — extrapolation to SWIR degrades results)
                    if wavelength < 1000.0:
                        R_atm = R_atm * path_corr_c * (wavelength / 550.0) ** path_corr_delta

                    # P3: NIR/SWIR transmittance correction from dark vegetation
                    # t_corr > 1 means retrieved > expected → T too low → increase T
                    if wavelength >= 700.0 and path_corr_t != 1.0:
                        t_corr_sqrt = np.sqrt(path_corr_t)
                        T_down = T_down * t_corr_sqrt
                        T_up = T_up * t_corr_sqrt

                    # Combined two-way transmittance for reliability checks
                    T_scat = T_down * T_up

                    # Check LUT reliability at this wavelength
                    scalar_lut = np.ndim(T_scat) == 0
                    band_reliable = True

                    if scalar_lut:
                        # Adaptive T_scat threshold: compare with expected
                        # scattering-only transmittance. If T_scat is much
                        # lower, gas absorption has contaminated the LUT.
                        wl_um = wavelength / 1000.0
                        tau_r = (0.008569 * wl_um**(-4) *
                                 (1 + 0.0113 * wl_um**(-2)) *
                                 pressure / 1013.25)
                        aod_scalar = (float(np.mean(p_aod_eff))
                                      if np.ndim(p_aod_eff) > 0
                                      else float(p_aod_eff))
                        tau_a = aod_scalar * (wavelength / 550.0)**(-1.3)
                        m_geom = (1.0 / np.cos(theta_s_rad) +
                                  1.0 / np.cos(np.radians(view_zenith)))
                        T_expected = np.exp(-(tau_r + tau_a) * 0.5 * m_geom)
                        # Wavelength-dependent factor: permissive at blue
                        # (gas-scattering coupling legitimately reduces T_scat)
                        # strict at NIR/SWIR (T_scat should match T_expected)
                        if wavelength < 550:
                            adapt_factor = 0.3
                        elif wavelength > 750:
                            adapt_factor = 0.6
                        else:
                            adapt_factor = 0.3 + 0.3 * (wavelength - 550) / 200
                        T_min_adaptive = adapt_factor * T_expected
                        effective_min = max(LUT_TSCAT_MIN, T_min_adaptive)

                        if T_scat < effective_min:
                            band_reliable = False
                            gs.verbose(
                                f"Band {band_num} ({wavelength:.1f} nm): "
                                f"LUT T_scat {T_scat:.4f} < "
                                f"{effective_min:.4f} (gas contamination)")
                        else:
                            # Spectral stability: if T_scat changes rapidly
                            # over ±10nm, we're at a gas band edge where LUT
                            # interpolation blends incompatible gas states
                            wl_lo = max(wavelength - 10.0,
                                        atm_lut.wavelengths[0])
                            wl_hi = min(wavelength + 10.0,
                                        atm_lut.wavelengths[-1])
                            _, Td_lo, Tu_lo, _ = atm_lut.interpolate(
                                wl_lo, p_aod_eff)
                            T_lo = Td_lo * Tu_lo
                            _, Td_hi, Tu_hi, _ = atm_lut.interpolate(
                                wl_hi, p_aod_eff)
                            T_hi = Td_hi * Tu_hi
                            T_max = max(T_lo, T_scat, T_hi)
                            T_min_v = max(min(T_lo, T_scat, T_hi), 1e-6)
                            if T_max / T_min_v > LUT_SPECTRAL_RATIO:
                                band_reliable = False
                                gs.verbose(
                                    f"Band {band_num} ({wavelength:.1f} nm): "
                                    f"T_scat unstable "
                                    f"({T_min_v:.3f}-{T_max:.3f})")

                    if not band_reliable:
                        refl_boa_band = np.full_like(refl_toa, np.nan)
                        sigma_rfl_band = np.full_like(refl_toa, np.nan)
                    else:
                        # 6S-style inversion: LUT already includes gas
                        # absorption at the per-pixel H2O (no tg_ratio needed)
                        T_total = T_down * T_up
                        y = (refl_toa - R_atm) / T_total
                        with np.errstate(divide='ignore', invalid='ignore'):
                            refl_boa_band = y / (1.0 + s * y)
                        # Clip valid values; NaN (nodata) stays NaN
                        finite = np.isfinite(refl_boa_band)
                        refl_boa_band[~finite] = np.nan
                        refl_boa_band[finite] = np.clip(
                            refl_boa_band[finite], -0.01, 1.5
                        )

                        # Mask unreliable LUT pixels (per-pixel AOD)
                        if not scalar_lut:
                            refl_boa_band[T_scat < LUT_TSCAT_MIN] = np.nan

                        # --- Improvement #2: In-loop adjacency correction ---
                        if adjacency_psf_km > 0 and pixel_size is not None:
                            try:
                                import adjacency as adj_module
                                r_env = adj_module.compute_environmental_reflectance(
                                    refl_boa_band, pixel_size, adjacency_psf_km)
                                T_dir = adj_module.compute_direct_transmittance(
                                    wavelength, p_aod, pressure, solar_zenith,
                                    view_zenith)
                                T_diff = np.clip(T_total - T_dir, 0.0, T_total)
                                correction = (T_diff * s * (refl_boa_band - r_env)
                                              / (1.0 - s * r_env + 1e-10))
                                refl_boa_band = refl_boa_band + correction
                                refl_boa_band = np.clip(
                                    np.nan_to_num(refl_boa_band, nan=np.nan),
                                    -0.01, 1.5)
                            except Exception as e:
                                if i == 0:
                                    gs.warning(f"In-loop adjacency failed: {e}")

                        # --- Improvement #4/#6: Per-band uncertainty ---
                        sigma_rfl_band = np.zeros_like(refl_boa_band)
                        if compute_uncertainty:
                            try:
                                import uncertainty as unc_module
                                sigma_rfl_band = unc_module.compute_reflectance_uncertainty(
                                    rad_toa_band, refl_boa_band,
                                    E0_band, d2, np.cos(theta_s_rad),
                                    T_down, T_up, s, R_atm,
                                    aod_sigma=0.04,
                                    atm_lut=atm_lut,
                                    wavelength=wavelength,
                                    aod=p_aod, h2o=p_h2o)
                            except Exception as e:
                                if i == 0:
                                    gs.warning(f"Uncertainty computation failed: {e}")

                # Store corrected band data in memory
                band_data_list.append((band, refl_boa_band))
                uncertainty_list.append((band, sigma_rfl_band))

                gs.percent(i, len(bands), 1)

            except Exception as e:
                gs.fatal(f"Error processing band {band_num}: {str(e)}")

        if not band_data_list:
            gs.fatal("No bands were successfully processed")

        # --- Improvement #3/#5: Surface prior regularization (MAP blend) ---
        if polish or compute_uncertainty:
            try:
                import surface_model as surface_model_module
                import uncertainty as unc_module

                wavelengths_for_prior = np.array(
                    [b['wavelength'] for b, _ in band_data_list])
                n_bp = len(band_data_list)
                ref_shape = band_data_list[0][1].shape

                # Stack reflectance into [n_bands, rows, cols]
                refl_stack = np.stack([d for _, d in band_data_list], axis=0)

                if refl_stack.ndim == 3:
                    surf_mdl = surface_model_module.SurfaceModel(wavelengths_for_prior)

                    # Build measurement variance from uncertainty + model discrepancy (#6)
                    sigma_obs2 = None
                    if compute_uncertainty and uncertainty_list:
                        unc_stack = np.stack([u for _, u in uncertainty_list], axis=0)
                        sigma_model = unc_module.compute_model_discrepancy(
                            wavelengths_for_prior)
                        # Add model discrepancy to measurement uncertainty
                        sigma_total_sq = (unc_stack ** 2 +
                                          sigma_model[:, np.newaxis, np.newaxis] ** 2)
                        sigma_obs2 = sigma_total_sq

                        # Update uncertainty_list with total including model discrepancy
                        sigma_total = np.sqrt(sigma_total_sq)
                        uncertainty_list = [
                            (binfo, sigma_total[idx])
                            for idx, (binfo, _) in enumerate(uncertainty_list)
                        ]

                    gs.message("Applying surface prior regularization (MAP)...")
                    refl_reg = surf_mdl.regularize(
                        refl_stack, sigma_obs2=sigma_obs2, weight=0.1)

                    # Write regularized data back to band_data_list
                    band_data_list = [
                        (binfo, refl_reg[idx])
                        for idx, (binfo, _) in enumerate(band_data_list)
                    ]
                    gs.message("  Surface prior regularization complete")

            except ImportError:
                gs.verbose("Surface model not available, skipping MAP regularization")
            except Exception as e:
                gs.warning(f"Surface prior regularization failed: {e}")

        # --- In-memory spectral polishing (before writing to rasters) ---
        if polish:
            import spectral_polish

            gs.message("Applying spectral polishing...")

            wavelengths_arr = np.array(
                [b['wavelength'] for b, _ in band_data_list])
            n_bands_p = len(band_data_list)

            # Compute per-band gas transmission for quality weighting
            tg_per_band = np.ones(n_bands_p)
            for idx, (binfo, _) in enumerate(band_data_list):
                try:
                    coefs_p = get_smac_parameters(
                        wavelength=binfo['wavelength'],
                        fwhm=binfo.get('fwhm', 10.0),
                        aerosol_model=aerosol_model,
                    )
                    tg_p = compute_gas_transmission(
                        coefs_p, solar_zenith, view_zenith,
                        pressure, water_vapor, ozone,
                    )
                    tg_per_band[idx] = float(tg_p)
                except Exception:
                    tg_per_band[idx] = 1.0

            quality_weights = spectral_polish.compute_quality_weights(
                wavelengths_arr, tg_per_band
            )

            # Stack band data: [n_bands, rows, cols]
            ref_shape = band_data_list[0][1].shape
            all_data = np.stack([d for _, d in band_data_list], axis=0)

            if all_data.ndim == 3:
                rows, cols = all_data.shape[1], all_data.shape[2]
                pixels = all_data.reshape(n_bands_p, -1).T
            else:
                # 1D bands (single row)
                pixels = all_data.T

            polished, pflags = spectral_polish.spectral_polish(
                pixels, wavelengths_arr,
                quality_weights=quality_weights,
                window=15, mad_threshold=2.5, replace=True,
            )

            n_flagged = np.sum(pflags)
            gs.message(f"Spectral polishing: {n_flagged} band-pixel values "
                       f"flagged and interpolated")

            # Write polished data back
            if all_data.ndim == 3:
                polished_3d = polished.T.reshape(n_bands_p, rows, cols)
            else:
                polished_3d = polished.T

            band_data_list = [
                (binfo, polished_3d[idx])
                for idx, (binfo, _) in enumerate(band_data_list)
            ]

        # Write corrected bands to GRASS rasters
        gs.message("Writing corrected bands...")
        for binfo, refl_data in band_data_list:
            band_num = binfo['band_num']
            wavelength = binfo['wavelength']
            fwhm = binfo.get('fwhm', 10.0)

            output_band = f"tmp_corr_{os.getpid()}_{band_num}"
            corrected_bands.append(output_band)

            ncols = refl_data.shape[1] if refl_data.ndim > 1 else refl_data.shape[0]
            with RasterRow(output_band, mode='w', mtype='DCELL',
                           overwrite=True) as r:
                for row_data in refl_data:
                    buf = Buffer((ncols,), mtype='DCELL')
                    buf[:] = row_data
                    r.put_row(buf)

            band_comment = f"Band {band_num}: {wavelength:.2f} nm"
            if fwhm:
                band_comment += f", FWHM: {fwhm:.2f} nm"
            gs.run_command('r.support', map=output_band,
                          title=band_comment, units="reflectance",
                          quiet=True)

        # Restore 3D region and combine bands
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        gs.message("Combining corrected bands into 3D raster...")
        gs.run_command('r.to.rast3', input=','.join(corrected_bands),
                      output=output_raster, overwrite=True)

        # --- Write uncertainty raster (Improvement #4) ---
        if compute_uncertainty and output_uncertainty and uncertainty_list:
            gs.message("Writing uncertainty bands...")
            unc_band_names = []
            for binfo, sigma_data in uncertainty_list:
                band_num = binfo['band_num']
                unc_name = f"tmp_unc_{os.getpid()}_{band_num}"
                unc_band_names.append(unc_name)

                ncols = (sigma_data.shape[1] if sigma_data.ndim > 1
                         else sigma_data.shape[0])
                with RasterRow(unc_name, mode='w', mtype='DCELL',
                               overwrite=True) as r:
                    for row_data in sigma_data:
                        buf = Buffer((ncols,), mtype='DCELL')
                        buf[:] = row_data
                        r.put_row(buf)

            gs.run_command('r.to.rast3', input=','.join(unc_band_names),
                          output=output_uncertainty, overwrite=True)
            gs.message(f"Uncertainty raster written: {output_uncertainty}")

            # Clean up uncertainty temp bands
            for name in unc_band_names:
                if gs.find_file(name, element='cell')['file']:
                    gs.run_command('g.remove', flags='f', type='raster',
                                  name=name, quiet=True)

        # Transfer metadata
        try:
            desc = ["Atmospheric Correction Metadata (libRadtran LUT):"]
            desc.append(f"Original raster: {input_raster}")
            desc.append("Method: libRadtran LUT (DISORT multiple-scattering)")
            desc.append(f"Solar Z: {solar_zenith}, View Z: {view_zenith}")
            desc.append(f"AOD: {aod}, Water Vapor: {water_vapor} g/cm2")
            desc.append(f"Ozone: {ozone} cm-atm, Pressure: {pressure} hPa")
            desc.append(f"Aerosol model: {aerosol_model}")
            desc.append("Measurement: Reflectance (Bottom of Atmosphere)")

            if wavelength_info:
                desc.append(f"Valid Bands: {len(corrected_bands)}")
                for j, band in enumerate(bands):
                    if j < len(corrected_bands):
                        wl_line = next(
                            (w for w in wavelength_info
                             if f"Band {band['band_num']}:" in w), None
                        )
                        if wl_line:
                            desc.append(wl_line)

            gs.run_command('r3.support', map=output_raster,
                          title=f"LUT corrected {input_raster}",
                          description="\n".join(desc),
                          source1="GRASS GIS i.hyper.smac module (libRadtran LUT)",
                          quiet=True)

            try:
                timestamp = gs.read_command('r3.timestamp',
                                            map=input_raster).strip()
                if timestamp:
                    gs.run_command('r3.timestamp', map=output_raster,
                                  date=timestamp)
            except Exception:
                pass

        except Exception as e:
            gs.warning(f"Could not transfer metadata: {str(e)}")

        gs.percent(1, 1, 1)
        gs.message(f"LUT atmospheric correction complete: {output_raster}")

    except Exception as e:
        gs.fatal(f"Error in LUT processing: {str(e)}")

    finally:
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        if not keep_temp:
            gs.message("Cleaning up temporary files...")
            for temp_band in temp_bands:
                if gs.find_file(temp_band, element='cell')['file']:
                    gs.run_command('g.remove', flags='f', type='raster',
                                  name=temp_band, quiet=True)
            for corr_band in corrected_bands:
                if gs.find_file(corr_band, element='cell')['file']:
                    gs.run_command('g.remove', flags='f', type='raster',
                                  name=corr_band, quiet=True)
        
        # Clean up OpenCL resources
        if accelerator:
            accelerator.cleanup()
            if gs.verbosity() > 1:
                gs.message("OpenCL resources cleaned up")


def _apply_spectral_polishing(output_raster, bands, solar_zenith, view_zenith,
                               pressure, water_vapor, ozone, aod, aerosol_model):
    """Apply spectral polishing to the corrected 3D raster in-place.

    Reads all corrected bands, runs moving-median outlier detection along
    the spectral axis, and overwrites flagged bands with interpolated values.
    """
    import spectral_polish

    gs.message("Applying spectral polishing...")

    wavelengths = np.array([b['wavelength'] for b in bands])
    n_bands = len(bands)

    # Compute per-band gas transmission for quality weighting
    tg_per_band = np.ones(n_bands)
    for idx, band in enumerate(bands):
        try:
            coefs = get_smac_parameters(
                wavelength=band['wavelength'],
                fwhm=band.get('fwhm', 10.0),
                aerosol_model=aerosol_model,
            )
            tg = compute_gas_transmission(
                coefs, solar_zenith, view_zenith,
                pressure, water_vapor, ozone,
            )
            tg_per_band[idx] = float(tg)
        except Exception:
            tg_per_band[idx] = 1.0

    quality_weights = spectral_polish.compute_quality_weights(
        wavelengths, tg_per_band
    )

    # Read all bands into memory [n_bands, rows, cols]
    all_data = []
    for band in bands:
        band_num = band['band_num']
        temp_name = f"tmp_polish_{os.getpid()}_{band_num}"
        gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
        gs.run_command('r3.to.rast', input=output_raster,
                       output=temp_name, overwrite=True, quiet=True)
        gs.run_command('g.rename',
                       raster=f"{temp_name}_00001,{temp_name}",
                       overwrite=True, quiet=True)
        with RasterRow(temp_name) as r:
            r.open('r')
            all_data.append(np.array(r))
        gs.run_command('g.remove', flags='f', type='raster',
                       name=temp_name, quiet=True)

    all_data = np.array(all_data)  # [n_bands, rows, cols]
    rows, cols = all_data.shape[1], all_data.shape[2]

    # Reshape to [n_pixels, n_bands] for spectral_polish
    pixels = all_data.reshape(n_bands, -1).T  # [n_pixels, n_bands]

    polished, flags = spectral_polish.spectral_polish(
        pixels, wavelengths, quality_weights=quality_weights,
        window=15, mad_threshold=2.5, replace=True,
    )

    n_flagged = np.sum(flags)
    gs.message(f"Spectral polishing: {n_flagged} band-pixel values flagged "
               f"and interpolated")

    # Reshape back and write corrected bands
    polished_3d = polished.T.reshape(n_bands, rows, cols)

    corrected_band_names = []
    for idx, band in enumerate(bands):
        band_num = band['band_num']
        out_name = f"tmp_polished_{os.getpid()}_{band_num}"
        corrected_band_names.append(out_name)
        ncols = cols
        with RasterRow(out_name, mode='w', mtype='DCELL', overwrite=True) as r:
            for row_data in polished_3d[idx]:
                buf = Buffer((ncols,), mtype='DCELL')
                buf[:] = row_data
                r.put_row(buf)

    # Re-pack into 3D raster
    gs.run_command('g.region', raster_3d=output_raster, quiet=True)
    band_list = ','.join(corrected_band_names)
    gs.run_command('r.to.rast3', input=band_list,
                   output=output_raster, overwrite=True, quiet=True)

    for name in corrected_band_names:
        gs.run_command('g.remove', flags='f', type='raster',
                       name=name, quiet=True)


def _apply_adjacency_correction(output_raster, bands, solar_zenith, view_zenith,
                                 aod, pressure, aerosol_model, psf_radius_km,
                                 force_regenerate=False):
    """Apply adjacency effect correction to the corrected 3D raster in-place.

    For each band, reads the BOA reflectance, computes the environmental
    reflectance via spatial averaging, and applies the Vermote et al. (1997)
    correction.
    """
    import adjacency
    import lut as lut_module

    gs.message("Applying adjacency effect correction...")

    # Get pixel size from GRASS region
    region = gs.region()
    pixel_size = (region['ewres'] + region['nsres']) / 2.0

    # Generate a LUT for scattering parameters
    phi = 0.0
    atm_lut = lut_module.AtmosphericLUT.get_or_generate(
        sza=solar_zenith, vza=view_zenith, phi=phi,
        pressure=pressure, aerosol_model=aerosol_model,
        force_regenerate=force_regenerate,
    )

    corrected_band_names = []
    for band in bands:
        band_num = band['band_num']
        wavelength = band['wavelength']

        if wavelength < 300 or wavelength > 4000:
            continue

        temp_name = f"tmp_adj_{os.getpid()}_{band_num}"
        gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
        gs.run_command('r3.to.rast', input=output_raster,
                       output=temp_name, overwrite=True, quiet=True)
        gs.run_command('g.rename',
                       raster=f"{temp_name}_00001,{temp_name}",
                       overwrite=True, quiet=True)

        with RasterRow(temp_name) as r:
            r.open('r')
            r_boa = np.array(r)

        # Get scattering parameters from LUT
        R_atm, T_down, T_up, s = atm_lut.interpolate(wavelength, aod)
        T_scat = T_down * T_up

        r_adj = adjacency.adjacency_correction(
            r_boa, T_scat, s, wavelength, aod, pressure,
            solar_zenith, view_zenith, pixel_size, psf_radius_km,
        )
        r_adj = np.clip(np.nan_to_num(r_adj, nan=np.nan), -0.01, 1.5)

        out_name = f"tmp_adj_out_{os.getpid()}_{band_num}"
        corrected_band_names.append(out_name)
        ncols = r_adj.shape[1] if r_adj.ndim > 1 else r_adj.shape[0]
        with RasterRow(out_name, mode='w', mtype='DCELL', overwrite=True) as r:
            for row_data in r_adj:
                buf = Buffer((ncols,), mtype='DCELL')
                buf[:] = row_data
                r.put_row(buf)

        gs.run_command('g.remove', flags='f', type='raster',
                       name=temp_name, quiet=True)

    # Re-pack into 3D raster
    gs.run_command('g.region', raster_3d=output_raster, quiet=True)
    band_list = ','.join(corrected_band_names)
    gs.run_command('r.to.rast3', input=band_list,
                   output=output_raster, overwrite=True, quiet=True)

    for name in corrected_band_names:
        gs.run_command('g.remove', flags='f', type='raster',
                       name=name, quiet=True)

    gs.message("Adjacency correction complete.")


def main():
    """Main function."""
    options, flags = gs.parser()
    
    input_raster = options['input']
    output_raster = options['output']
    dem = options['dem']
    keep_temp = flags['k']
    clear_lut_cache = flags['c']
    method = 'lut'
    
    # Get viewing geometry
    if options['solar_zenith']:
        solar_zenith = float(options['solar_zenith'])
    else:
        gs.fatal("Solar Zenith Angle is required. Please provide the solar_zenith parameter.")
    
    if options['solar_azimuth']:
        solar_azimuth = float(options['solar_azimuth'])
    else:
        gs.fatal("Solar Azimuth Angle is required. Please provide the solar_azimuth parameter.")
    
    if options['view_zenith']:
        view_zenith = float(options['view_zenith'])
    else:
        gs.fatal("View Zenith Angle is required. Please provide the view_zenith parameter.")
    
    if options['view_azimuth']:
        view_azimuth = float(options['view_azimuth'])
    else:
        gs.fatal("View Azimuth Angle is required. Please provide the view_azimuth parameter.")
    
    # Get atmospheric parameters
    if options['pressure']:
        pressure = float(options['pressure'])
    else:
        pressure = estimate_pressure_from_dem(dem)
    
    # Initialize default AOD value and map
    aod = 0.15  # Typical clear atmosphere
    aod_map = None  # Initialize aod_map to None
    retrieved_alpha = None  # Angstrom exponent from Dark Target

    if options['aod'] and options['aod'].strip():  # Check for non-empty string
        try:
            aod = float(options['aod'])
            gs.message(f"Using provided aod value: {aod}")
        except ValueError:
            gs.message("AOD not provided, estimating from hyperspectral data...")
    else:
        # Estimate AOD if not provided
        aod_map, aod, retrieved_alpha = estimate_aod(
            input_raster=input_raster,
            dem=dem,
            method='auto',
            solar_zenith=solar_zenith,
            view_zenith=view_zenith,
            solar_azimuth=solar_azimuth,
            view_azimuth=view_azimuth,
            pressure=pressure,
            verbose=gs.verbosity() > 1
        )
        gs.message(f"Estimated AOD @ 550nm: {aod:.3f}")
        if retrieved_alpha is not None:
            gs.message(f"Retrieved Angstrom exponent: {retrieved_alpha:.2f}")

    # Fix 3: Clamp AOD to upper bound
    if aod > AOD_MAX:
        gs.warning(f"AOD value {aod:.3f} exceeds maximum {AOD_MAX}. Clamping to {AOD_MAX}.")
        aod = AOD_MAX

    # Estimate ozone if not provided
    ozone = 0.4  # Typical value
    ozone_map = None  # Initialize ozone_map to None
    
    if options['ozone'] and options['ozone'].strip():
        try:
            ozone = float(options['ozone'])
            gs.message(f"Using provided ozone value: {ozone} cm-atm")
        except:
            gs.message("Ozone not provided, estimating from hyperspectral data...")    
    else:
        # Estimate ozone using Chappuis band method
        ozone_map, ozone_du = o3.estimate_ozone(
            input_raster=input_raster,
            method='chappuis',
            verbose=gs.verbosity() > 1
        )
        # Convert from DU to cm-atm (1 DU = 0.001 cm-atm)
        ozone = ozone_du * 0.001
        gs.message(f"Estimated total column ozone: {ozone_du:.1f} DU ({ozone:.3f} cm-atm)")
                        
    # Log the ozone value being used
    gs.message(f"Ozone: {ozone:.3f} cm-atm")

    # Initialize default water vapor content.
    gs.message("WVC: Estimating water vapor content...")
    water_vapor = 2.0  # g/cm² - typical mid-latitude value
    wvc_map = None  # Initialize wvc_map to None
    wvc_uncertainty = 0.1  # Default uncertainty
    wvc_method = 'joint'  # Default method
    
    # Parse water_vapor parameter - can be numeric value or method name
    if options['water_vapor']:
        wv_input = options['water_vapor'].strip().lower()
        if wv_input in ['joint', '940nm', '1130nm', 'average']:
            # It's a method specification
            wvc_method = wv_input
            gs.message(f"Using WVC estimation method: {wvc_method}")
        else:
            # It's a numeric water vapor value
            try:
                water_vapor = float(wv_input)
                gs.message(f"Using provided water vapor value: {water_vapor} g/cm²")
            except ValueError:
                gs.warning(f"Invalid water_vapor value '{wv_input}', using default estimation")
    
    if wvc_method != 'joint' and not options['water_vapor']:
        # If method wasn't specified and we're not using a numeric value, estimate WVC
        gs.message("Water vapor not provided, estimating from hyperspectral data...")
        try:
            # Estimate WVC from hyperspectral data using enhanced joint retrieval
            wvc_map, water_vapor, wvc_uncertainty = estimate_wvc(
                input_raster=input_raster,
                dem=dem,
                method=wvc_method,  # Use specified method
                solar_zenith=solar_zenith,
                view_zenith=view_zenith,
                verbose=gs.verbosity() > 0
            )
            gs.message(f"Estimated water vapor content: {water_vapor:.2f} g/cm² (±{wvc_uncertainty:.2f})")
        except Exception as e:
            gs.warning(f"Failed to estimate water vapor from data: {str(e)}")
            gs.warning("Falling back to default water vapor value")
    
    # Get band information
    gs.message(f"Processing {input_raster}...")
    bands = get_all_band_wavelengths(input_raster)
    
    gs.message(f"Found {len(bands)} bands")
    gs.message(f"Wavelength range: {bands[0]['wavelength']:.1f} - {bands[-1]['wavelength']:.1f} nm")
    
    # Print atmospheric parameters
    gs.message("=" * 60)
    gs.message("Atmospheric Parameters:")
    gs.message(f"  Method: lut (libRadtran DISORT, enhanced gas absorption)")
    if options['water_vapor'] and options['water_vapor'].strip().lower() in ['joint', '940nm', '1130nm', 'average']:
        gs.message(f"  WVC Retrieval: {wvc_method} (optimal estimation)")
    else:
        gs.message(f"  WVC Retrieval: estimation from data")
    gs.message(f"  Solar zenith: {solar_zenith:.1f}°")
    gs.message(f"  AOD @ 550nm: {aod:.3f}")
    gs.message(f"  Water vapor: {water_vapor:.2f} g/cm²")
    if 'wvc_uncertainty' in locals():
        gs.message(f"  WVC uncertainty: ±{wvc_uncertainty:.2f} g/cm²")
    gs.message(f"  Ozone: {ozone:.2f} cm-atm")
    gs.message(f"  Pressure: {pressure:.1f} hPa")

    aerosol_model = options.get('aerosol_model', 'continental')
    angstrom_alpha = float(options['angstrom']) if options.get('angstrom') else retrieved_alpha
    apply_polish = flags['p']
    adjacency_psf = float(options.get('adjacency_psf', 0))
    do_uncertainty = flags.get('u', False)
    unc_output = options.get('output_uncertainty', '')

    gs.message(f"  Aerosol model: {aerosol_model}")
    if options.get('angstrom'):
        gs.message(f"  Angstrom exponent: {angstrom_alpha:.2f} (user override)")
    elif angstrom_alpha is not None:
        gs.message(f"  Angstrom exponent: {angstrom_alpha:.2f} (retrieved from Dark Target)")
    else:
        gs.message(f"  Angstrom exponent: native Shettle (default for {aerosol_model})")

    if apply_polish:
        gs.message("  Spectral polishing: ENABLED")
    if adjacency_psf > 0:
        gs.message(f"  Adjacency correction: PSF radius = {adjacency_psf:.1f} km")
    if do_uncertainty:
        gs.message("  Uncertainty estimation: ENABLED")
        if unc_output:
            gs.message(f"  Uncertainty output: {unc_output}")

    gs.message("=" * 60)

    # Apply LUT correction (only method)
    apply_lut_correction(
        input_raster, output_raster, bands,
        aod, water_vapor, ozone, pressure,
        solar_zenith, solar_azimuth,
        view_zenith, view_azimuth,
        aerosol_model, keep_temp,
        aod_map=aod_map,
        wvc_map=wvc_map,
        ozone_map=ozone_map,
        dem=dem,
        polish=apply_polish,
        angstrom_alpha=angstrom_alpha,
        force_regenerate=clear_lut_cache,
        adjacency_psf_km=adjacency_psf,
        compute_uncertainty=do_uncertainty,
        output_uncertainty=unc_output if unc_output else None,
        opencl_device=options.get('opencl_device', 'auto'),
        opencl_memory=int(options.get('opencl_memory', 1024)),
        parallel_lut=options.get('parallel_lut', 'auto'),
        smart_lut=options.get('smart_lut', 'auto')
    )

    return 0

if __name__ == "__main__":
    sys.exit(main())
