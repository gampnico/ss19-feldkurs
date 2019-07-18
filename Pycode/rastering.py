import rasterio
import rasterio.plot
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import georaster
from osgeo import osr
import scipy.io
import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import matplotlib.pylab as plt
import re
import os
import glob


#
# mat = scipy.io.loadmat('./rasters/q3m.mat')
# ref = osr.SpatialReference()
# # Import from Proj.4
#
# # Or EPSG code
# ref.ImportFromEPSG(31254)
# tly_max = 237155.63 + 2486
# tly_min = 237155.63 - 2846
# tlx_min = 79029.64 - 2486
# tlx_max = 79029.64 + 2486
#
# # trans_q3 = (tlx, 10, 0, tly, 0, -10)
#
# array = mat["q3m"]["fclim"][0][0] * 10 ** 9
# array[array < 1] = 0
# # My image array
# lat = mat["q3m"]["yn"][0][0]
# lon = mat["q3m"]["xn"][0][0]
# # For each pixel I know it's latitude and longitude.
# # As you'll see below you only really need the coordinates of
# # one corner, and the resolution of the file.
#
# xmin, ymin, xmax, ymax = [tlx_min, tly_min, tlx_max, tly_max]
# nrows, ncols = np.shape(array)
# xres = (xmax - xmin) / float(ncols)
# yres = (ymax - ymin) / float(nrows)
#
# geotransform = (xmin, xres, 0, ymax, 0, -yres)
# print(geotransform)
# # That's (top left x, w-e pixel resolution, rotation (0 if North is up),
# #         top left y, rotation (0 if North is up), n-s pixel resolution)
# # I don't know why rotation is in twice???
#
# output_raster = gdal.GetDriverByName('GTiff').Create('myraster.tif', ncols,
#                                                      nrows, 1,
#                                                      gdal.GDT_Float64)  #
# output_raster.GetRasterBand(1).WriteArray(
#     array)  # Writes my array to the raster
# output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
# srs = osr.SpatialReference()  # Establish its coordinate encoding
# srs.ImportFromEPSG(31254)  # This one specifies WGS84 lat long.
# output_raster.SetProjection(
#     srs.ExportToWkt())  # Exports the coordinate system to the file
# output_raster.FlushCache()
#
# drv = gdal.GetDriverByName('GTiff')
# ds_in = gdal.Open('./rasters/q3r.asc')
# ds_out = drv.CreateCopy('./rasters/q3r.tif', ds_in)
# srs = osr.SpatialReference()
# srs.ImportFromEPSG(31254)
# ds_out.SetProjection(srs.ExportToWkt())
# ds_in = None
# ds_out = None
#
# gtif = gdal.Open("./rasters/q3r.tif")
# band = gtif.GetRasterBand(1)
# bandArray = band.ReadAsArray()
# print(bandArray)
# print(np.max(bandArray))
# proj4 = '+proj=stere ...'
# im = georaster.SingleBandRaster.from_array(array, geotransform, proj4)
# im.save_geotiff('my_filename.tif')

def asc_convert(file_name):
    """
    Args:
        file_name (str): name of text file to be converted to ASCII raster.
        Can use search
    Returns:


    """
    txt_in = "./rasters/asc_loose/" + str(file_name) + ".txt"
    asc_out = "./rasters/asc_loose/" + str(file_name) + ".asc"
    f = open(txt_in, 'r')
    filedata = f.read()
    f.close()
    newdata = re.sub('\s+', ' ', filedata).strip()
    newdata = re.sub(',', ' ', newdata).strip()
    newdata = re.sub(';', '\n', newdata).strip()
    # newdata = filedata.replace(r" +", r" ")
    print(newdata)
    f = open(asc_out, 'w')
    f.write(newdata)
    f.close()


# import rasterio
# from rasterio.merge import merge
# from rasterio.plot import show
# import glob
# import os
#
# dirpath = r"./rasters"
#
# out_fp = r"./rasters/mosaic.tif"
#
# search_criteria = "q*.tif"
# q = os.path.join(dirpath, search_criteria)
# print(q)
# dem_fps = glob.glob(q)
# print(dem_fps)
#
# src_files_to_mosaic = []
#
# for fp in dem_fps:
#     src = rasterio.open(fp)
#     src_files_to_mosaic.append(src)
#
# mosaic, out_trans = merge(src_files_to_mosaic)
# out_meta = src.meta.copy()
#
# out_meta.update({
#     "driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2],
#     "transform": out_trans, "crs": "+proj=utm +zone=35 +ellps=GRS80
#     +units=m "
#                                    "+no_defs "})
# with rasterio.open(out_fp, "w", **out_meta) as dest:
#     dest.write(mosaic)
# mosaic, out_trans = gdal.Warp(src_files_to_mosaic)

def asc_to_tif(q_num):
    """
    Converts ASCII raster file into tiff image format
    Args:
        q_num (str): filename of fclim raster. Can use search_for_files()
            to generate list of valid files.
    """
    drv = gdal.GetDriverByName('GTiff')
    in_path = "./rasters/asc_loose/q" + str(q_num) + ".asc"
    out_path = "./rasters/gen_tif/q" + str(q_num) + ".tif"
    ds_in = gdal.Open(in_path)
    ds_out = drv.CreateCopy(out_path, ds_in)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(31254)
    ds_out.SetProjection(srs.ExportToWkt())
    ds_in = None
    ds_out = None


# q_num_list = ["125", "250", "375", "500", "625", "750", "875"]
# for num in q_num_list:
#     asc_to_tif(num)
def read_file(filename):
    """
    Authors:
        rutgerhofste (2015)

    Args:

    Returns:

    """
    filehandle = gdal.Open(filename)
    band1 = filehandle.GetRasterBand(1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    Z = band1.ReadAsArray()
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize, ysize, geotransform, geoproj, Z


def write_file(filename, geotransform, geoprojection, data):
    """
    Authors:
        rutgerhofste (2015)

    Args:

    Returns:

    """
    (x, y) = data.shape
    file_format = "GTiff"
    driver = gdal.GetDriverByName(file_format)
    # you can change the data format but be sure to be able to store negative
    # values including -9999
    dst_datatype = gdal.GDT_Float64
    dst_ds = driver.Create(filename, y, x, 1, dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds = None
    return 1


def search_for_files(search_criteria="q*.tif"):
    dirpath = r"./rasters"
    # search_criteria = "*_new.tif"
    q = os.path.join(dirpath, search_criteria)
    print(q)
    dem_fps = glob.glob(q)
    return dem_fps


def tiff_transform(q_num_list, weight_list):
    """
    Authors:
        rutgerhofste (2015), Nicolas Gampierakis (2019)

    Args:

    Returns:

    """
    firstrun = 1
    sum_weight = sum(weight_list)
    for q_num in q_num_list:
        in_path = "./rasters/gen_tif/q" + str(q_num) + ".tif"
        out_path = "./rasters/proc_tif/q" + str(q_num) + ".tif"
        index = q_num_list.index(q_num)
        weight = weight_list[index]
        if firstrun == 1:
            [xsize, ysize, geotransform, geoproj, Z] = read_file(in_path)

            Z = Z * weight / sum_weight

            Z[np.isnan(Z)] = -9999
            new_geotransform = (
                geotransform[0], geotransform[1], geotransform[2],
                geotransform[3], geotransform[4], -1 * geotransform[5])
            print(new_geotransform)
            # set the new GeoTransform, effectively flipping the image
            write_file(out_path, new_geotransform, geoproj, Z)
