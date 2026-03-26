/*
 * r3d_slice.c — Fast Z-slice extraction from a GRASS 3D raster using
 *               Rast3d_get_block() in no-cache (tile-bulk) mode.
 *
 * Compiled to lib/libr3dslice.so and called from lib/r3d_slice.py.
 *
 * Function:
 *   int r3d_read_z_slice(const char *name3d, const char *mapset3d,
 *                        int z, int ncols, int nrows, double *buf);
 *
 * Parameters:
 *   name3d    - 3D raster map name
 *   mapset3d  - mapset (NULL or "" → auto-find)
 *   z         - depth index (0 = bottom = first band = shortest wavelength)
 *   ncols     - expected number of columns (must match map region)
 *   nrows     - expected number of rows    (must match map region)
 *   buf       - caller-allocated double[nrows * ncols], row-major
 *
 * Returns 0 on success, -1 on error.
 *
 * Block layout from Rast3d_get_block():
 *   buf[z_off * ny * nx + y * nx + x]  with z_off=0, ny=nrows, nx=ncols
 * i.e. buf[y * ncols + x] after the call — standard row-major 2-D layout.
 */

#include <stdio.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster3d.h>

/* Query the spatial dimensions of a 3D raster without opening it.
 * Returns 0 on success; out params are set to the map's region.
 */
int r3d_map_dims(const char *name3d, const char *mapset3d,
                 int *ncols_out, int *nrows_out, int *depths_out)
{
    RASTER3D_Region reg;
    const char *ms = (mapset3d && mapset3d[0]) ? mapset3d
                                                : G_find_raster3d(name3d, "");
    if (!ms) {
        fprintf(stderr, "r3d_map_dims: map '%s' not found\n", name3d);
        return -1;
    }
    Rast3d_init_defaults();
    if (!Rast3d_read_region_map(name3d, ms, &reg)) {
        fprintf(stderr, "r3d_map_dims: Rast3d_read_region_map failed\n");
        return -1;
    }
    *ncols_out  = reg.cols;
    *nrows_out  = reg.rows;
    *depths_out = reg.depths;
    return 0;
}

int r3d_read_z_slice(const char *name3d, const char *mapset3d,
                     int z, int ncols, int nrows, double *buf)
{
    RASTER3D_Region map_region;
    RASTER3D_Map   *map;
    const char     *ms;

    /* Resolve mapset */
    if (mapset3d && mapset3d[0] != '\0') {
        ms = mapset3d;
    } else {
        ms = G_find_raster3d(name3d, "");
        if (!ms) {
            fprintf(stderr, "r3d_read_z_slice: map '%s' not found\n", name3d);
            return -1;
        }
    }

    /* Read the map's native region so we can open the full cube */
    if (!Rast3d_read_region_map(name3d, ms, &map_region)) {
        fprintf(stderr, "r3d_read_z_slice: Rast3d_read_region_map failed for '%s'\n", name3d);
        return -1;
    }

    /* Validate dimensions */
    if (map_region.cols != ncols || map_region.rows != nrows) {
        fprintf(stderr,
                "r3d_read_z_slice: dimension mismatch: map is %d×%d, caller expects %d×%d\n",
                map_region.cols, map_region.rows, ncols, nrows);
        return -1;
    }
    if (z < 0 || z >= map_region.depths) {
        fprintf(stderr,
                "r3d_read_z_slice: z=%d out of range [0, %d)\n",
                z, map_region.depths);
        return -1;
    }

    /* Open with full native region + no-cache for bulk tile reads */
    Rast3d_init_defaults();
    map = Rast3d_open_cell_old(name3d, ms, &map_region,
                               RASTER3D_TILE_SAME_AS_FILE,
                               RASTER3D_NO_CACHE);
    if (!map) {
        fprintf(stderr, "r3d_read_z_slice: Rast3d_open_cell_old failed for '%s'\n", name3d);
        return -1;
    }

    /*
     * Rast3d_get_block(map, x0, y0, z0, nx, ny, nz, buf, type)
     * Reads a contiguous block starting at (x0,y0,z0) of size (nx,ny,nz).
     * With nz=1 we get a single Z-slice.
     * Block layout: buf[(z_off)*ny*nx + y*nx + x]
     * → with z_off=0: buf[y*ncols + x] — standard row-major.
     */
    Rast3d_get_block(map, 0, 0, z, ncols, nrows, 1, buf, DCELL_TYPE);

    Rast3d_close(map);
    return 0;
}
