import numpy as np
from osgeo import gdal, osr

# Input and output files
input_file = "band_111.raw"
output_file = "band_111.tif"

# Read the raw data
data = np.fromfile(input_file, dtype='float32')
data = data.reshape((732, 607))  # Reshape to 2D array

# Create output GeoTIFF
driver = gdal.GetDriverByName('GTiff')
ds = driver.Create(output_file, 607, 732, 1, gdal.GDT_Float32)

# Set the geotransform
# (xmin, xres, 0, ymax, 0, -yres)
xmin = 406800
ymax = 2918190
xres = (431430 - 406800) / 607
yres = (2918190 - 2890950) / 732
ds.SetGeoTransform((xmin, xres, 0, ymax, 0, -yres))

# Set the projection (UTM zone 44N, WGS84)
srs = osr.SpatialReference()
srs.ImportFromEPSG(32644)  # WGS 84 / UTM zone 44N
ds.SetProjection(srs.ExportToWkt())

# Write the data
ds.GetRasterBand(1).WriteArray(data)
ds.GetRasterBand(1).SetNoDataValue(-9999)
ds = None  # Close the file
