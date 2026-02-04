#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate SMAC coefficients for hyperspectral bands.

This script generates SMAC coefficient files for a range of wavelengths
suitable for hyperspectral atmospheric correction.

Usage:
    python generate_hyperspectral_coefs.py --start 400 --end 2500 --step 10
    python generate_hyperspectral_coefs.py --sensor PRISMA
    python generate_hyperspectral_coefs.py --wavelengths 450,550,650,850

Author: i.hyper.smac GRASS GIS module
"""

import os
import sys
import argparse
import numpy as np

# Add lib directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from smac_coef_generator import generate_smac_coefficients, SMACCoefficients


# Predefined sensor band configurations
SENSOR_BANDS = {
    'PRISMA': {
        'description': 'PRISMA Hyperspectral (ASI)',
        'wavelengths': np.arange(400, 2505, 10),  # 400-2500nm, ~10nm bands
        'fwhm': 10.0,
    },
    'AVIRIS': {
        'description': 'AVIRIS Classic (NASA)',
        'wavelengths': np.arange(400, 2500, 10),
        'fwhm': 10.0,
    },
    'AVIRIS_NG': {
        'description': 'AVIRIS Next Generation (NASA)',
        'wavelengths': np.arange(380, 2510, 5),  # 5nm resolution
        'fwhm': 5.0,
    },
    'HYPERION': {
        'description': 'EO-1 Hyperion',
        'wavelengths': np.concatenate([
            np.arange(356, 1058, 10),  # VNIR
            np.arange(852, 2577, 10),  # SWIR
        ]),
        'fwhm': 10.0,
    },
    'ENMAP': {
        'description': 'EnMAP (DLR)',
        'wavelengths': np.arange(420, 2450, 6.5),
        'fwhm': 6.5,
    },
    'EMIT': {
        'description': 'EMIT (NASA/JPL)',
        'wavelengths': np.arange(381, 2493, 7.4),
        'fwhm': 7.4,
    },
    'TANAGER': {
        'description': 'Planet Tanager',
        'wavelengths': np.arange(400, 2500, 5),
        'fwhm': 5.0,
    },
    'GENERIC_VNIR': {
        'description': 'Generic VNIR Hyperspectral (400-1000nm)',
        'wavelengths': np.arange(400, 1001, 5),
        'fwhm': 5.0,
    },
    'GENERIC_FULL': {
        'description': 'Generic Full Range (400-2500nm)',
        'wavelengths': np.arange(400, 2501, 10),
        'fwhm': 10.0,
    },
}

# Aerosol models
AEROSOL_MODELS = ['continental', 'maritime', 'urban', 'desert']


def generate_coefficients_for_wavelengths(wavelengths, fwhm, aerosol_type,
                                          output_dir, verbose=False):
    """
    Generate SMAC coefficients for a list of wavelengths.

    Args:
        wavelengths: Array of wavelengths in nm
        fwhm: Full width at half maximum in nm
        aerosol_type: Aerosol model type
        output_dir: Output directory for coefficient files
        verbose: Enable verbose output

    Returns:
        List of generated coefficient file paths
    """
    os.makedirs(output_dir, exist_ok=True)

    generated_files = []
    total = len(wavelengths)

    print(f"\nGenerating {total} coefficient files...")
    print(f"Aerosol model: {aerosol_type}")
    print(f"Output directory: {output_dir}\n")

    for i, wl in enumerate(wavelengths):
        filename = f"coef_{wl:.0f}nm_{aerosol_type.upper()}.dat"
        filepath = os.path.join(output_dir, filename)

        print(f"[{i+1}/{total}] Generating {filename}...", end=' ', flush=True)

        try:
            coef = generate_smac_coefficients(
                wavelength_nm=wl,
                fwhm_nm=fwhm,
                aerosol_type=aerosol_type,
                output_file=filepath,
                verbose=False  # Suppress per-band verbose output
            )
            generated_files.append(filepath)
            print("OK")

        except Exception as e:
            print(f"FAILED: {e}")

            # Create a fallback coefficient file with analytical values
            if verbose:
                print(f"  Creating fallback coefficients...")

            coef = create_analytical_coefficients(wl, aerosol_type)
            coef.to_file(filepath)
            generated_files.append(filepath)

    return generated_files


def create_analytical_coefficients(wavelength_nm: float,
                                   aerosol_type: str = 'continental') -> SMACCoefficients:
    """
    Create analytical SMAC coefficients when libRadtran is not available.

    These are approximate coefficients based on physical models.

    Args:
        wavelength_nm: Wavelength in nanometers
        aerosol_type: Aerosol model type

    Returns:
        SMACCoefficients with analytical values
    """
    coef = SMACCoefficients()

    wl_um = wavelength_nm / 1000.0

    # Rayleigh optical thickness (Hansen & Travis 1974)
    coef.taur = 0.008569 * wl_um**(-4) * (1 + 0.0113 * wl_um**(-2))
    coef.sr = coef.taur

    # Aerosol optical depth scaling (Ångström law, α=1.3)
    coef.a0taup = 0.0
    coef.a1taup = (wavelength_nm / 550.0) ** (-1.3)

    # Aerosol properties by type
    aerosol_props = {
        'continental': (0.89, 0.65),
        'maritime': (0.98, 0.72),
        'urban': (0.82, 0.62),
        'desert': (0.92, 0.72),
    }
    coef.wo, coef.gc = aerosol_props.get(aerosol_type, (0.89, 0.65))

    # Spherical albedo (simplified)
    coef.a0s = coef.taur * 0.5
    coef.a1s = 0.2
    coef.a2s = -0.05
    coef.a3s = 0.0

    # Transmission (simplified)
    coef.a0T = 1.0
    coef.a1T = -0.15
    coef.a2T = 0.0
    coef.a3T = -0.1

    # O3 absorption (Chappuis band: 450-700nm)
    if 450 < wavelength_nm < 700:
        # Peak around 600nm
        o3_strength = np.exp(-((wavelength_nm - 600) / 80) ** 2)
        coef.ao3 = -0.1 * o3_strength
        coef.no3 = 1.0
    else:
        coef.ao3 = 0.0
        coef.no3 = 1.0

    # H2O absorption
    # Strong bands at 720, 820, 940, 1130, 1380, 1900 nm
    h2o_bands = [720, 820, 940, 1130, 1380, 1900]
    h2o_widths = [20, 30, 50, 50, 60, 100]
    h2o_strengths = [0.01, 0.02, 0.1, 0.15, 0.3, 0.5]

    coef.ah2o = 0.0
    for band, width, strength in zip(h2o_bands, h2o_widths, h2o_strengths):
        if abs(wavelength_nm - band) < width * 2:
            contrib = strength * np.exp(-((wavelength_nm - band) / width) ** 2)
            coef.ah2o -= contrib

    coef.nh2o = 0.5

    # Phase function (Henyey-Greenstein approximation)
    g = coef.gc
    coef.a0P = (1 - g**2) / (1 + g)**1.5
    coef.a1P = 0.0
    coef.a2P = 0.0
    coef.a3P = 0.0
    coef.a4P = 0.0

    return coef


def main():
    parser = argparse.ArgumentParser(
        description='Generate SMAC coefficients for hyperspectral bands',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate for a specific sensor
  python generate_hyperspectral_coefs.py --sensor PRISMA

  # Generate for a wavelength range
  python generate_hyperspectral_coefs.py --start 400 --end 1000 --step 10

  # Generate for specific wavelengths
  python generate_hyperspectral_coefs.py --wavelengths 450,550,650,850

  # Generate for all aerosol models
  python generate_hyperspectral_coefs.py --sensor AVIRIS --all-aerosols

Available sensors:
""" + '\n'.join([f"  {k}: {v['description']}" for k, v in SENSOR_BANDS.items()])
    )

    parser.add_argument('--sensor', type=str, choices=list(SENSOR_BANDS.keys()),
                        help='Predefined sensor configuration')
    parser.add_argument('--start', type=float, help='Start wavelength (nm)')
    parser.add_argument('--end', type=float, help='End wavelength (nm)')
    parser.add_argument('--step', type=float, default=10.0,
                        help='Wavelength step (nm)')
    parser.add_argument('--wavelengths', type=str,
                        help='Comma-separated list of wavelengths (nm)')
    parser.add_argument('--fwhm', type=float, default=10.0,
                        help='Full width at half maximum (nm)')
    parser.add_argument('--aerosol', type=str, default='continental',
                        choices=AEROSOL_MODELS,
                        help='Aerosol model type')
    parser.add_argument('--all-aerosols', action='store_true',
                        help='Generate for all aerosol models')
    parser.add_argument('--output', type=str, default='../COEFS',
                        help='Output directory')
    parser.add_argument('--analytical', action='store_true',
                        help='Use analytical formulas (no libRadtran)')
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose output')

    args = parser.parse_args()

    # Determine wavelengths
    if args.sensor:
        config = SENSOR_BANDS[args.sensor]
        wavelengths = config['wavelengths']
        fwhm = config['fwhm']
        print(f"\nSensor: {args.sensor}")
        print(f"Description: {config['description']}")
    elif args.wavelengths:
        wavelengths = np.array([float(w) for w in args.wavelengths.split(',')])
        fwhm = args.fwhm
    elif args.start and args.end:
        wavelengths = np.arange(args.start, args.end + args.step, args.step)
        fwhm = args.fwhm
    else:
        parser.error("Specify --sensor, --wavelengths, or --start/--end")

    print(f"\nWavelength range: {wavelengths.min():.0f} - {wavelengths.max():.0f} nm")
    print(f"Number of bands: {len(wavelengths)}")
    print(f"FWHM: {fwhm:.1f} nm")

    # Determine aerosol models
    if args.all_aerosols:
        aerosol_types = AEROSOL_MODELS
    else:
        aerosol_types = [args.aerosol]

    # Resolve output directory
    output_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), args.output)
    )

    # Generate coefficients
    all_files = []
    for aerosol in aerosol_types:
        if args.analytical:
            # Use analytical formulas
            print(f"\nGenerating analytical coefficients for {aerosol}...")
            subdir = os.path.join(output_dir, aerosol.upper())
            os.makedirs(subdir, exist_ok=True)

            for wl in wavelengths:
                filename = f"coef_{wl:.0f}nm_{aerosol.upper()}.dat"
                filepath = os.path.join(subdir, filename)
                coef = create_analytical_coefficients(wl, aerosol)
                coef.to_file(filepath)
                all_files.append(filepath)

            print(f"Generated {len(wavelengths)} files in {subdir}")
        else:
            # Use libRadtran
            subdir = os.path.join(output_dir, aerosol.upper())
            files = generate_coefficients_for_wavelengths(
                wavelengths, fwhm, aerosol, subdir, args.verbose
            )
            all_files.extend(files)

    print(f"\n{'='*60}")
    print(f"Generation complete!")
    print(f"Total files generated: {len(all_files)}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
