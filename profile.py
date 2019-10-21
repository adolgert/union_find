# coding: utf-8
import sys;
import numpy as np
import unittest
import gdal
from gdalconst import *
import raster_stats as raster_stats
reload(raster_stats)


def test_file():
    filename='34418039.tif'
    ds=gdal.Open( filename.encode("utf8"), GA_ReadOnly )
    ds_array = ds.ReadAsArray()
    return raster_stats.find_clusters_time(1,ds_array)
