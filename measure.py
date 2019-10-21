import sys
import gdal

from gdalconst import *


if __name__ == '__main__':
    if len(sys.argv)>1:
        filename=sys.argv[1]
    else:
        filename='34418039.tif'
    print filename, type(filename)
    dataset = gdal.Open( filename.encode("utf8"), GA_ReadOnly )
    ds_numpy_array = dataset.ReadAsArray()
