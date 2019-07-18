# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 21:33:20 2015

@author: rutgerhofste

This script will replace values in a raster and save the result in geotiff
format

"""

from osgeo import gdal
import numpy as np

firstrun = 1


def readFile(filename):
    filehandle = gdal.Open(filename)
    band1 = filehandle.GetRasterBand(1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    Z = band1.ReadAsArray()
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize, ysize, geotransform, geoproj, Z


def writeFile(filename, geotransform, geoprojection, data):
    (x, y) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    # you can change the dataformat but be sure to be able to store negative
    # values including -9999
    dst_datatype = gdal.GDT_Float32
    dst_ds = driver.Create(filename, y, x, 1, dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    return 1


#
# pathname = "./rasters/qmid.tif"
# writefilename = "./rasters/qmid_new.tif"
#
# if firstrun == 1:
#     [xsize, ysize, geotransform, geoproj, Z] = readFile(pathname)
#
# # Set large negative values to -9999
# Z[Z < 0] = -9999
# # Set very small values to -9999
# Z[Z < 1*10**-18] = -9999
# # Choose your preference: (comment either rule)
# # Z[Z == -9999] = np.nan
# # Or
# Z[np.isnan(Z)] = -9999
#
# writeFile(writefilename, geotransform, geoproj, Z)
import math

x = (78900 - 77450)
y = (236140 - 235700)
a = math.sqrt(x ** 2 + y ** 2)

print(x)
print(y)
print(a)

