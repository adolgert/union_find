# coding: utf-8
import sys;
from datetime import datetime
import numpy as np
import unittest
import logging
from argparse import ArgumentParser
from glob import glob
from itertools import permutations
from scipy.stats import stats
import yaml # module py27-yaml
import gdal
from gdalconst import *
import raster_stats as raster_stats
reload(raster_stats)


logger=logging.getLogger('timing')



def test_file():
    filename='34418039.tif'
    ds=gdal.Open( filename.encode("utf8"), GA_ReadOnly )
    ds_array = ds.ReadAsArray()
    return ds_array



def resize_array(arr,N):
    '''
    We want to return an array of a new size that may be larger
    or smaller than the original array.
    '''
    assert(len(arr.shape)==2)

    nrr = np.zeros(N*N, dtype=arr.dtype).reshape((N,N))
    # Loop over the number of copies of arr into new array
    for i in range((N-1)/arr.shape[0] + 1):
        irange=[i*arr.shape[0],(i+1)*arr.shape[0]]
        if irange[1]>N:
            irange[1]=N

        for j in range((N-1)/arr.shape[1] + 1):
            jrange=[j*arr.shape[1],(j+1)*arr.shape[1]]
            if jrange[1]>N:
                jrange[1]=N
            logger.debug('tile %d,%d: %s %s, %d %d' % (i,j,
                         str(irange), str(jrange),
                         irange[1]-irange[0], jrange[1]-jrange[0]))
            nrr[irange[0]:irange[1],jrange[0]:jrange[1]] = \
                arr[0:irange[1]-irange[0],0:jrange[1]-jrange[0]]
    return nrr



class TestResize(unittest.TestCase):
    def test_identity(self):
        arr=np.array([1,2,3,1,2,3,4,5,6]).reshape((3,3))
        r=resize_array(arr,3)
        logging.debug('identity: %s %s' % (str(arr), str(r)))
        self.assertTrue((arr==r).all())

    def test_smaller(self):
        arr=np.array([1,2,3,1,2,3,4,5,6]).reshape((3,3))
        r=resize_array(arr,2)
        self.assertTrue((r==np.array([1,2,1,2]).reshape((2,2))).all())

    def test_twice(self):
        arr=np.array([1,2,3,1,2,3,4,5,6]).reshape((3,3))
        r=resize_array(arr,6)
        self.assertTrue((arr==r[0:3,0:3]).all())
        self.assertTrue((arr==r[3:6,0:3]).all())
        self.assertTrue((arr==r[0:3,3:6]).all())
        self.assertTrue((arr==r[3:6,3:6]).all())

    def test_larger(self):
        arr=np.array([1,2,3,1,2,3,4,5,6]).reshape((3,3))
        r=resize_array(arr,4)
        self.assertTrue((r[0:3,0:3]==arr).all())
        self.assertTrue((r[3:4,0:3]==arr[0:1,0:3]).all())
        self.assertTrue((r[0:3,3:4]==arr[0:3,0:1]).all())
        self.assertEqual(r[3,3],arr[0,0])





def run_tests(tests,N=500,run_cnt=30):
    arr=test_file()
    marr=resize_array(arr,N)

    all_results=dict()
    for t in tests:
        logger.info('Timing %s for %d runs' % (t,run_cnt))
        times=list()
        g=getattr(raster_stats,t)
        f=lambda x: g(1,x)
        for i in range(run_cnt):
            times.append(f(marr))

        all_results[t]=times
        results={ 'size' : (N,N), 'test' : t, 'times' : times, 'when' : datetime.now() }
        with open('%s.yaml' % t,'w') as f:
            yaml.dump(results,f)
            logger.info('wrote results to %s' % f.name)

        logger.info('avg %gs' % (np.average(times)/10**9))
        if run_cnt>1:
            logger.info('stdev %gs' % (np.std(times)/10**9))

    del marr

    if run_cnt>1 and len(tests)>1:
        for a,b in permutations(tests,2):
            print a,b
            print '\t', stats.ttest_ind(all_results[a],all_results[b])


def heap_tests(tests,N=500,run_cnt=1):
    arr=test_file()
    marr=resize_array(arr,N)

    all_results=dict()
    for t in tests:
        logger.info('Heap test %s for %d runs' % (t,run_cnt))
        times=list()
        g=getattr(raster_stats,t)

        for i in range(run_cnt):
            times.append(g(marr))

    del marr
    #raster_stats.dump_heap()



def show_tests():
    test_files=glob('*.yaml')
    results=dict()
    for filename in test_files:
        with open(filename,'r') as f:
            blob=yaml.load(f)
            results[blob['test']]=blob

    for name,vals in results.items():
        times=vals['times']
        print name, vals['size'], vals['when']
        print 'avg', np.average(times)/10**9
        print 'stdev', np.std(times)/10**9


def suite():
    suite=unittest.TestSuite()
    loader=unittest.TestLoader()
    suite.addTests(loader.loadTestsFromTestCase(TestResize))
    return suite



if __name__ == '__main__':
    parser = ArgumentParser(description='profile union-find algorithm')
    parser.add_argument('--size', dest='size', type=int, default=500,
                        help='The dimensions of the NxN array')
    parser.add_argument('--reps', dest='reps', type=int, default=1,
                        help='The number of times to repeat each run.')
    parser.add_argument('--heap', dest='heap', action='store_true',
                        help='just run without timing')
    parser.add_argument('--func', dest='func', type=str, default='all',
                        help='which function to test, or all')
    parser.add_argument('--log', dest='log', type=str, default='INFO',
               help='Logging level: TRACE, DEBUG, INFO, WARN, ERROR')
    parser.add_argument('--test', dest='test', action='store_true',
               help='Whether to run unit tests on this file.')
    args = parser.parse_args()

    allowed_debug=['TRACE', 'DEBUG', 'INFO', 'WARN', 'ERROR']
    if not args.log in allowed_debug:
        print 'Debugging level should be one of %s' % (str(allowed_debug),)

    logging.basicConfig(level=getattr(logging,args.log))

    if args.test:
        unittest.TextTestRunner(verbosity=2).run(suite())
        sys.exit(0)

    if args.func == 'all':
        tests=['find_clusters_time','find_clusters_pointer_time','find_clusters_remap_time']
    else:
        if not args.func in dir(raster_stats):
            logging.error('Could not find %s in raster_stats module.' % \
                              args.func)
            sys.exit(1)
        tests=[args.func]

    if args.heap:
        heap_tests(tests,args.size,args.reps)
    else:
        run_tests(tests,args.size,args.reps)

    
