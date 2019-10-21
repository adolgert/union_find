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


def unique():
    filename='34418039.tif'
    ds=gdal.Open( filename.encode("utf8"), GA_ReadOnly )
    ds_array = ds.ReadAsArray()
    return raster_stats.unique_values(ds_array)
    

class KnownArrays(unittest.TestCase):
    def test_one(self):
        t=np.array([1,1,1,1],dtype=np.uint8).reshape((2,2))
        r=raster_stats.find_clusters(t)
        self.assertEqual(len(r),1)
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0][0],0)
        self.assertEqual(results[0][1],1)
        self.assertEqual(results[0][2],2)
        self.assertEqual(results[0][3],3)

    def test_two(self):
        t=np.array([1,1,2,2],dtype=np.uint8).reshape((2,2))
        r=raster_stats.find_clusters(t)
        self.assertEqual(len(r),2)
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0][0],0)
        self.assertEqual(results[0][1],1)
        self.assertEqual(results[1][0],2)
        self.assertEqual(results[1][1],3)

    def test_three(self):
        t=np.array([1,1,2,6],dtype=np.uint8).reshape((2,2))
        r=raster_stats.find_clusters(t)
        self.assertEqual(len(r),3)
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0][0],0)
        self.assertEqual(results[0][1],1)
        self.assertEqual(results[1][0],2)
        self.assertEqual(results[2][0],3)

    def test_bigger(self):
        t=np.array([1,2,3, 1,2,3, 3,2,1],dtype=np.uint8).reshape((3,3))
        r=raster_stats.find_clusters(t)
        self.assertEqual(len(r),5)
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0],[0,3])
        self.assertEqual(results[1],[1,4,7])
        self.assertEqual(results[2],[2,5])
        self.assertEqual(results[3],[6])
        self.assertEqual(results[4],[8])
        
        
    def test_subarray(self):
        '''
        The raster_stats library copies a pointer from the numpy
        array. In subarrays, that pointer still points to the 
        larger array.
        '''
        t=np.array([1,2,3, 1,2,3, 3,2,1],dtype=np.uint8).reshape((3,3))
        # Will a smaller subset of this array be passed as a contiguous
        # c pointer? Should be [[1,2],[1,2]] not [[1,2],[3,1]]
        r=raster_stats.find_clusters(t[0:2,0:2])
        self.assertEqual(len(r),2,"The wrapper misreads subarrays of numpy arrays.")
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0],[0,2])
        self.assertEqual(results[1],[1,3])


    def test_subarray_copy(self):
        '''
        Because Python can't understand subarrays of numpy arrays, copy
        the subarray to a new array.
        '''
        t=np.array([1,2,3, 1,2,3, 3,2,1],dtype=np.uint8).reshape((3,3))
        # Will a smaller subset of this array be passed as a contiguous
        # c pointer? Should be [[1,2],[1,2]] not [[1,2],[3,1]]
        r=raster_stats.find_clusters(t[0:2,0:2].copy())
        self.assertEqual(len(r),2,"The wrapper misreads subarrays of numpy arrays.")
        results=list()
        for x in r:
            results.append(x)
        self.assertEqual(results[0],[0,2])
        self.assertEqual(results[1],[1,3])
        
        


def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(KnownArrays)
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=4).run(suite())
    