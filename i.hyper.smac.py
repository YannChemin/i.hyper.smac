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
# % type: double
# % required: no
# % description: Water vapor content (g/cm²)
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
# % key: method
# % type: string
# % required: no
# % options: simple,libradtran,lut
# % answer: simple
# % description: Atmospheric correction method to use (simple, libradtran, or lut)
# % guisection: Method
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

# %flag
# % key: k
# % description: Keep temporary bands
# % guisection: Optional
# %end

# %flag
# % key: g
# % description: Generate SMAC coefficients from libRadtran (requires libRadtran and scipy)
# % guisection: Advanced
# %end

import sys
import os
import numpy as np
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


AOD_MAX = 1.5  # Upper bound for AOD validation
TRANSMISSION_THRESHOLD = 0.10  # Minimum two-way transmission


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
                                   ozone_map=None, dem=None):
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
                # Approximate E0 using 5778K blackbody scaled to solar constant
                # F(λ) = Ω_sun × B(λ, T_sun), result in mW/(m² nm)
                wl_m = wavelength * 1e-9  # nm to meters
                hc_kT = 6.626e-34 * 2.998e8 / (wl_m * 1.381e-23 * 5778.0)
                B = 2 * 6.626e-34 * (2.998e8)**2 / (wl_m**5 * (np.exp(hc_kT) - 1.0))
                E0_band = 6.794e-5 * B * 1e-6  # W/m³ → mW/(m² nm)
            
            # Convert radiance to TOA reflectance
            # Formula: refl_toa = (π * L_toa * d²) / (E0 * cos(θs))
            refl_toa_band = (np.pi * rad_toa_band * d2) / (E0_band * np.cos(theta_s_rad))
            
            # Keep negative values as-is for debugging (indicates calibration issues)
            # refl_toa_band = np.clip(refl_toa_band, 0.0, 1.0)
            
            try:
                # Get SMAC coefficients from libRadtran
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

                # Check two-way transmission — mask opaque bands
                band_T = compute_band_transmission(
                    coefs, solar_zenith, view_zenith,
                    p_pres, p_aod, p_wvc, p_o3
                )

                if np.ndim(band_T) == 0:
                    # Scalar path (no maps)
                    if band_T < TRANSMISSION_THRESHOLD:
                        gs.verbose(
                            f"Band {band_num} ({wavelength:.1f} nm): "
                            f"transmission {band_T:.3f} < {TRANSMISSION_THRESHOLD}, masking as NaN"
                        )
                        refl_boa_band = np.full_like(refl_toa_band, np.nan)
                    else:
                        aod_band = compute_blue_aod_taper(wavelength, p_aod)
                        refl_boa_band = smac.smac_inv(
                            r_toa=refl_toa_band,
                            tetas=solar_zenith,
                            phis=solar_azimuth,
                            tetav=view_zenith,
                            phiv=view_azimuth,
                            pressure=p_pres,
                            taup550=aod_band,
                            uo3=p_o3,
                            uh2o=p_wvc,
                            coef=coefs
                        )
                        refl_boa_band = np.clip(refl_boa_band, -0.01, 1.5)
                else:
                    # Per-pixel path
                    opaque_mask = band_T < TRANSMISSION_THRESHOLD
                    if np.all(opaque_mask):
                        gs.verbose(
                            f"Band {band_num} ({wavelength:.1f} nm): "
                            f"all pixels opaque, masking as NaN"
                        )
                        refl_boa_band = np.full_like(refl_toa_band, np.nan)
                    else:
                        aod_band = compute_blue_aod_taper(wavelength, p_aod)
                        refl_boa_band = smac.smac_inv(
                            r_toa=refl_toa_band,
                            tetas=solar_zenith,
                            phis=solar_azimuth,
                            tetav=view_zenith,
                            phiv=view_azimuth,
                            pressure=p_pres,
                            taup550=aod_band,
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
                         ozone_map=None, dem=None):
    """Apply atmospheric correction using libRadtran LUT.

    Uses full multiple-scattering from libRadtran DISORT via a precomputed
    look-up table for R_atm, T_scat, and s.  Gas absorption is handled
    separately via SMAC coefficients (supports per-pixel WVC).

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
    """
    import lut as lut_module

    gs.message("Applying libRadtran LUT atmospheric correction...")

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

    # Generate or load cached LUT
    gs.message("Generating/loading atmospheric LUT...")
    atm_lut = lut_module.AtmosphericLUT.get_or_generate(
        sza=solar_zenith, vza=view_zenith, phi=phi,
        pressure=pressure, aerosol_model=aerosol_model,
        wl_min=wl_min, wl_max=wl_max, wl_step=10,
    )

    # Get wavelength information from input raster
    input_info = gs.read_command('r3.info', flags='h', map=input_raster)
    wavelength_info = [line.strip() for line in input_info.split('\n')
                       if 'Band' in line and 'nm' in line]

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

    theta_s_rad = np.radians(solar_zenith)

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

            # Select per-pixel or scalar atmospheric parameters
            p_pres = pressure_array if pressure_array is not None else pressure
            p_aod = aod_array if aod_array is not None else aod
            p_wvc = wvc_array if wvc_array is not None else water_vapor
            p_o3 = ozone_array if ozone_array is not None else ozone

            try:
                # Get SMAC coefficients for gas transmission only
                coefs = get_smac_parameters(
                    wavelength=wavelength, fwhm=fwhm,
                    sza=solar_zenith, vza=view_zenith,
                    aod_550=aod, water_vapor=water_vapor,
                    ozone=ozone, pressure=pressure,
                    aerosol_model=aerosol_model,
                    verbose=gs.verbosity() > 0
                )

                # Compute gas-only transmission
                tg = compute_gas_transmission(
                    coefs, solar_zenith, view_zenith, p_pres, p_wvc, p_o3
                )

                # Mask opaque bands
                if np.ndim(tg) == 0:
                    if tg < TRANSMISSION_THRESHOLD:
                        gs.verbose(
                            f"Band {band_num} ({wavelength:.1f} nm): "
                            f"gas transmission {tg:.3f} < {TRANSMISSION_THRESHOLD}, "
                            f"masking as NaN"
                        )
                        refl_boa_band = np.full_like(refl_toa, np.nan)
                    else:
                        # Get scattering params from LUT
                        R_atm, T_scat, s = atm_lut.interpolate(wavelength, p_aod)

                        # Inversion formula
                        numerator = refl_toa - R_atm * tg
                        denominator = tg * T_scat + numerator * s
                        with np.errstate(divide='ignore', invalid='ignore'):
                            refl_boa_band = numerator / denominator
                        refl_boa_band = np.clip(
                            np.nan_to_num(refl_boa_band, nan=0.0), -0.01, 1.5
                        )
                else:
                    # Per-pixel path
                    opaque_mask = tg < TRANSMISSION_THRESHOLD
                    if np.all(opaque_mask):
                        gs.verbose(
                            f"Band {band_num} ({wavelength:.1f} nm): "
                            f"all pixels opaque, masking as NaN"
                        )
                        refl_boa_band = np.full_like(refl_toa, np.nan)
                    else:
                        R_atm, T_scat, s = atm_lut.interpolate(wavelength, p_aod)

                        numerator = refl_toa - R_atm * tg
                        denominator = tg * T_scat + numerator * s
                        with np.errstate(divide='ignore', invalid='ignore'):
                            refl_boa_band = numerator / denominator
                        refl_boa_band = np.clip(
                            np.nan_to_num(refl_boa_band, nan=0.0), -0.01, 1.5
                        )
                        refl_boa_band[opaque_mask] = np.nan

                # Write corrected band
                output_band = f"tmp_corr_{os.getpid()}_{band_num}"
                corrected_bands.append(output_band)

                ncols = refl_boa_band.shape[1] if refl_boa_band.ndim > 1 else refl_boa_band.shape[0]
                with RasterRow(output_band, mode='w', mtype='DCELL', overwrite=True) as r:
                    for row_idx, row_data in enumerate(refl_boa_band):
                        buf = Buffer((ncols,), mtype='DCELL')
                        buf[:] = row_data
                        r.put_row(buf)

                band_comment = f"Band {band_num}: {wavelength:.2f} nm"
                if fwhm:
                    band_comment += f", FWHM: {fwhm:.2f} nm"
                gs.run_command('r.support', map=output_band,
                              title=band_comment, units="reflectance",
                              quiet=True)

                gs.percent(i, len(bands), 1)

            except Exception as e:
                gs.fatal(f"Error processing band {band_num}: {str(e)}")

        if not corrected_bands:
            gs.fatal("No bands were successfully processed")

        # Restore 3D region and combine bands
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        gs.message("Combining corrected bands into 3D raster...")
        gs.run_command('r.to.rast3', input=','.join(corrected_bands),
                      output=output_raster, overwrite=True)

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


def main():
    """Main function."""
    options, flags = gs.parser()
    
    input_raster = options['input']
    output_raster = options['output']
    dem = options['dem']
    keep_temp = flags['k']
    generate_coefs = flags['g']
    method = options.get('method', 'simple')
    
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

    if options['aod'] and options['aod'].strip():  # Check for non-empty string
        try:
            aod = float(options['aod'])
            gs.message(f"Using provided aod value: {aod}")
        except ValueError:
            gs.message("AOD not provided, estimating from hyperspectral data...")
    else:
        # Estimate AOD if not provided
        aod_map, aod = estimate_aod(
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

    if options['water_vapor']:
        water_vapor = float(options['water_vapor'])
    else:
        gs.message("Water vapor not provided, estimating from hyperspectral data...")
        try:
            # Estimate WVC from hyperspectral data
            wvc_map, water_vapor = estimate_wvc(
                input_raster=input_raster,
                dem=dem,
                method='average',
                solar_zenith=solar_zenith,
                view_zenith=view_zenith,
                verbose=gs.verbosity() > 0
            )
            gs.message(f"Estimated water vapor content: {water_vapor:.2f} g/cm²")
                            
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
    gs.message(f"  Method: {method}")
    gs.message(f"  Solar zenith: {solar_zenith:.1f}°")
    gs.message(f"  AOD @ 550nm: {aod:.3f}")
    gs.message(f"  Water vapor: {water_vapor:.2f} g/cm²")
    gs.message(f"  Ozone: {ozone:.2f} cm-atm")
    gs.message(f"  Pressure: {pressure:.1f} hPa")
    
    if generate_coefs and method not in ('libradtran', 'lut'):
        gs.warning("The -g flag (generate coefficients from libRadtran) is only used "
                   "with method=libradtran. Ignoring -g flag.")
        generate_coefs = False

    aerosol_model = options.get('aerosol_model', 'continental')

    if method == 'libradtran':
        sensor_type = options.get('sensor', '').upper()
        visibility = float(options['visibility']) if options.get('visibility') else None

        gs.message(f"  Sensor: {sensor_type}")
        gs.message(f"  Aerosol model: {aerosol_model}")
        if visibility:
            gs.message(f"  Visibility: {visibility} km")
        if generate_coefs:
            gs.message("  Coefficient generation: ENABLED (libRadtran + scipy fitting)")

    if method == 'lut':
        gs.message(f"  Aerosol model: {aerosol_model}")

    gs.message("=" * 60)
    
    # Apply the selected correction method
    if method == 'simple':
        apply_smac_correction_simple(
            input_raster, output_raster, bands,
            aod, water_vapor, ozone, pressure,
            solar_zenith, solar_azimuth,
            view_zenith, view_azimuth,
            keep_temp
        )
    elif method == 'libradtran':
        apply_smac_correction_libradtran(
            input_raster, output_raster, bands,
            aod, water_vapor, ozone, pressure,
            solar_zenith, solar_azimuth,
            view_zenith, view_azimuth,
            sensor_type, visibility,
            aerosol_model, keep_temp,
            generate_coefs,
            aod_map=aod_map,
            wvc_map=wvc_map,
            ozone_map=ozone_map,
            dem=dem,
        )
    elif method == 'lut':
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
        )
    else:
        gs.fatal(f"Unknown method: {method}. Choose 'simple', 'libradtran', or 'lut'.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
