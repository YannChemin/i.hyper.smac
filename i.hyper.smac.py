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
# % key: aerosol_type
# % type: string
# % required: no
# % options: continental,maritime,urban,desert
# % answer: continental
# % description: Aerosol type for atmospheric correction
# % guisection: Atmospheric
# %end

# %option
# % key: method
# % type: string
# % required: no
# % options: simple,libradtran
# % answer: simple
# % description: Atmospheric correction method to use (simple or libradtran)
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

import sys
import os
import tempfile
import numpy as np
import grass.script as gs
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
            
            # Apply the correction
            expr = f"{temp_band_corr} = ({temp_band} - {corr['rho_atm']:.6f}) / {corr['t_total']:.6f}"
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
            description = (
                "Atmospherically corrected using SMAC method\n"
                f"Original band: {band_wl_info if band_wl_info is not None else 'N/A'}\n"
                f"Solar Z: {solar_zenith}°, View Z: {view_zenith}°\n"
                f"AOD: {aod}, Water Vapor: {water_vapor} g/cm², Ozone: {ozone} cm-atm\n"
                f"Wavelength: {band_wavelength} nm"
            )
            if 'fwhm' in band:
                description += f", FWHM: {band['fwhm']} nm"
            # Use r.support to add metadata to the specific band
            gs.run_command('r.support',
                          map=temp_band_corr,
                          title=f"SMAC corrected {input_raster} band {band_num}",
                          description=description,
                          source1="GRASS GIS i.hyper.smac module",
                          history=timestamp,
                          semantic_label=f"band_{band_num}",
                          vdatum="WGS84",
                          quiet=True)
            
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
                                   sensor_type, aot_lut=None, visibility=None,
                                   aerosol_type='continental', keep_temp=False):
    """Apply libradtran-based SMAC atmospheric correction.
    
    Args:
        input_raster (str): Input 3D raster name
        output_raster (str): Output 3D raster name
        bands (list): List of band information dictionaries
        aod (float): Aerosol Optical Depth at 550nm
        water_vapor (float): Water vapor content (g/cm²)
        ozone (float): Ozone content (cm-atm)
        pressure (float): Atmospheric pressure (hPa)
        solar_zenith (float): Solar zenith angle (degrees)
        solar_azimuth (float): Solar azimuth angle (degrees)
        view_zenith (float): View zenith angle (degrees)
        view_azimuth (float): View azimuth angle (degrees)
        sensor_type (str): Sensor type (e.g., 'AVIRIS', 'PRISMA')
        aot_lut (str, optional): Path to AOT look-up table
        visibility (float, optional): Visibility in km
        aerosol_type (str): Type of aerosol model
        keep_temp (bool): Whether to keep temporary files
    """
    
    gs.message("Applying libradtran-based SMAC atmospheric correction...")
    temp_dir = Path(tempfile.mkdtemp(prefix='smac_libradtran_'))
    
    try:
        # Process each band
        output_bands = []
        
        for band in bands:
            band_num = band['band_num']
            wavelength = band['wavelength']
            
            # Skip bands outside the valid range for libradtran
            if wavelength < 300 or wavelength > 4000:
                gs.warning(f"Skipping band {band_num} with wavelength {wavelength} nm (outside 300-4000 nm range)")
                continue
                
            gs.message(f"Processing band {band_num}: {wavelength} nm")
            
            # Extract single band
            # Extract band from 3D raster
            input_band = f"tmp_input_{os.getpid()}_{band_num}"
            
            # Set the 3D region to the specific band (using band_num + 0.1 to ensure top > bottom)
            gs.run_command('g.region', t=band_num + 0.1, b=band_num, quiet=True)
            
            # Extract the band using the mask
            gs.run_command('r3.mapcalc',
                          expression=f"{input_band} = {input_raster}",
                          overwrite=True,
                          quiet=True)
            
            # Convert the 3D raster to 2D with overwrite
            gs.run_command('r3.to.rast',
                         input=input_band,
                         output=input_band,
                         overwrite=True,
                         quiet=True)
            
            # The output will be named output_map_00001
            input_file = f"{input_band}_00001"
            
            # Rename the output file to the desired name
            gs.run_command('g.rename',
                         raster=f"{input_file},{input_band}",
                         overwrite=True,
                         quiet=True)
            
            try:
                # Get SMAC parameters from libRadtran
                smac_params = get_smac_parameters(
                    wavelength=wavelength,
                    fwhm=band.get('fwhm', 10.0),  # Default to 10nm if not specified
                    sza=solar_zenith,
                    aod_550=aod,
                    water_vapor=water_vapor,
                    ozone=ozone,  # This will be the estimated or provided ozone value
                    surface_albedo=0.1,  # Initial guess, will be updated
                    aerosol_type=aerosol_type,
                    verbose=gs.verbosity() > 0
                )
                
                # Apply atmospheric correction using SMAC parameters
                # rho_surface = (L_toa - L_p) / (T_dir * T_dif) - s * rho_surface
                # Solving for rho_surface: rho_surface = (L_toa - L_p) / (T_dir * T_dif + s * (L_toa - L_p))
                
                corrected_band = f"{input_band}_corr"
                expr = (
                    f"{corrected_band} = float({input_band} - {smac_params['path_radiance']}) / "
                    f"({smac_params['direct_transmittance']} * {smac_params['diffuse_transmittance']} + "
                    f"{smac_params['spherical_albedo']} * ({input_band} - {smac_params['path_radiance']}))"
                )
                
                gs.run_command('r.mapcalc', expression=expr, overwrite=True, quiet=True)
                output_bands.append(corrected_band)
                
            except Exception as e:
                gs.fatal(f"Error processing band {band_num}: {str(e)}")
        
        # Restore 3D g.region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        # Combine corrected bands back into a 3D raster
        gs.message("Combining corrected bands...")
        gs.run_command('r.to.rast3', input=','.join(output_bands), output=output_raster, overwrite=True)
            
    except Exception as e:
        gs.fatal(f"Error in libradtran processing: {str(e)}")
        
    finally:
        # Restore 3D g.region
        gs.run_command('g.region', raster_3d=input_raster, quiet=True)

        # Clean up temporary files
        if not keep_temp:
            gs.message("Cleaning up temporary files...")
            temp_bands = gs.parse_command('g.list', type='raster', 
                                         pattern=f"{output_raster}_band*", 
                                         mapset='.')
            if temp_bands:
                gs.run_command('g.remove', flags='f', type='raster', 
                              name=','.join(temp_bands), quiet=True)
            
            temp_bands_corr = gs.parse_command('g.list', type='raster', 
                                             pattern=f"{output_raster}_band*_corr", 
                                             mapset='.')
            if temp_bands_corr:
                gs.run_command('g.remove', flags='f', type='raster', 
                              name=','.join(temp_bands_corr), quiet=True)
    
    gs.message(f"Libradtran-based atmospheric correction complete: {output_raster}")

def main():
    """Main function."""    
    options, flags = gs.parser()
    
    input_raster = options['input']
    output_raster = options['output']
    dem = options['dem']
    keep_temp = flags['k']
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
            verbose=gs.verbosity() > 1
        )
        gs.message(f"Estimated AOD @ 550nm: {aod:.3f}")
    
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
    
    if method == 'libradtran':
        sensor_type = options.get('sensor', '').upper()
        visibility = float(options['visibility']) if options.get('visibility') else None
        aerosol_type = options.get('aerosol_type', 'continental') # continental is default
        
        gs.message(f"  Sensor: {sensor_type}")
        gs.message(f"  Aerosol type: {aerosol_type}")
        if visibility:
            gs.message(f"  Visibility: {visibility} km")
    
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
            aerosol_type, keep_temp
        )
    else:
        gs.fatal(f"Unknown method: {method}. Choose 'simple' or 'libradtran'.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
