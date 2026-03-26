"""
r3d_slice.py — ctypes wrapper for libr3dslice.so

Provides read_z_slice() to bulk-read a single Z-depth from a GRASS 3D raster
without spawning any subprocesses (no g.region / r3.to.rast / g.rename).

Usage::

    from lib.r3d_slice import read_z_slice

    # z = band_num - 1  (z=0 is bottom = first band = shortest wavelength)
    band_data = read_z_slice("my_r3_map", "", z=0, ncols=800, nrows=600)
    # Returns float64 ndarray shape (nrows, ncols), NaN where GRASS null.
"""

import ctypes
import os
import numpy as np

# ---------------------------------------------------------------------------
# Load the shared library
# ---------------------------------------------------------------------------
_so_path = os.path.join(os.path.dirname(__file__), "libr3dslice.so")
try:
    _lib = ctypes.CDLL(_so_path)
except OSError as exc:
    raise ImportError(
        f"Cannot load {_so_path}. Run 'make' in lib/ first.\n{exc}"
    ) from exc

# int r3d_map_dims(const char *name3d, const char *mapset3d,
#                  int *ncols_out, int *nrows_out, int *depths_out)
_lib.r3d_map_dims.restype  = ctypes.c_int
_lib.r3d_map_dims.argtypes = [
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]

# int r3d_read_z_slice(const char *name3d, const char *mapset3d,
#                      int z, int ncols, int nrows, double *buf)
_lib.r3d_read_z_slice.restype  = ctypes.c_int
_lib.r3d_read_z_slice.argtypes = [
    ctypes.c_char_p,          # name3d
    ctypes.c_char_p,          # mapset3d  (NULL → auto-find)
    ctypes.c_int,             # z
    ctypes.c_int,             # ncols
    ctypes.c_int,             # nrows
    ctypes.POINTER(ctypes.c_double),  # buf  (double[nrows*ncols])
]

_dims_cache: dict = {}


def map_dims(name3d: str, mapset3d: str = "") -> tuple:
    """Return (ncols, nrows, depths) for a GRASS 3D raster (cached)."""
    key = (name3d, mapset3d)
    if key not in _dims_cache:
        nc, nr, nd = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        ret = _lib.r3d_map_dims(
            name3d.encode(),
            mapset3d.encode() if mapset3d else b"",
            ctypes.byref(nc), ctypes.byref(nr), ctypes.byref(nd),
        )
        if ret != 0:
            raise RuntimeError(f"r3d_map_dims failed for '{name3d}'")
        _dims_cache[key] = (nc.value, nr.value, nd.value)
    return _dims_cache[key]


def read_z_slice(name3d: str, mapset3d: str, z: int,
                 ncols: int = 0, nrows: int = 0) -> np.ndarray:
    """Read one Z-depth from a GRASS 3D raster into a 2-D float64 array.

    Parameters
    ----------
    name3d   : GRASS 3D raster map name
    mapset3d : mapset name, or "" to auto-search
    z        : depth index — 0 = bottom = first band (shortest wavelength)
    ncols    : columns (0 = auto from map dims)
    nrows    : rows    (0 = auto from map dims)

    Returns
    -------
    np.ndarray, shape (nrows, ncols), dtype float64
        GRASS null cells are represented as NaN.
    """
    if ncols == 0 or nrows == 0:
        _nc, _nr, _ = map_dims(name3d, mapset3d)
        if ncols == 0:
            ncols = _nc
        if nrows == 0:
            nrows = _nr
    buf = np.empty(nrows * ncols, dtype=np.float64)
    ptr = buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    ms_bytes = mapset3d.encode() if mapset3d else b""

    ret = _lib.r3d_read_z_slice(
        name3d.encode(),
        ms_bytes,
        ctypes.c_int(z),
        ctypes.c_int(ncols),
        ctypes.c_int(nrows),
        ptr,
    )
    if ret != 0:
        raise RuntimeError(
            f"r3d_read_z_slice failed (code {ret}) for map '{name3d}' z={z}"
        )

    return buf.reshape(nrows, ncols)
