# Raster Stats Union Find

In this directory are two sets of files, one for Daniel's fractal
measurements and one for a C++ implementation of union-find, called
from Python.

Daniel's files: Fractals.py, RSSR.py, ClusterScript.py, GISDelaunay.py
Union-find:
  timing.py - Runs versions of union-find in order to see what's faster.
  traster.py - Unit tests on union-find. Examples of use.
  io_ppm.{h,cpp} - Writes PPM files, as a double-check to see if data is correct.
  io_geotiff.{h,cpp} - Reads geotiff files from C++.
  raster_times.{h,cpp} - Uses Boost.Chrono to time a function.
  raster_wrap.cpp - The Boost.Python wrapper on the C++ functions.
  timing_harness.{h,cpp} - A class to help Python time C++ functions.
  main.cpp - Pure C++ program, just for figuring out code.
  cluster.{h,cpp} - Union-find
    - find_clusters_pair, assumes an input multi-array, uses pair(i,j).
    - find_clusters, uses a single size_t for each (i,j) entry.
    - find_clusters_twopass, combines stages of algorithm into one loop.
    - find_clusters_pointer, returns a shared_ptr to the results.
    - find_clusters_remap, same cluster construction, but builds result map better.
  cluster_tbb.cpp - Union-find with TBB
    - clusters_tbb0, each (i,j) is a set of size_t. Parallel algorithm.
  cluster_generic.hpp - Union-find with generic templates and TBB

Requirements:

python2.7 - Would work with others, but this is the default.
numpy and scipy
py27-gdal is used to read geotiff in the script.
libgeotiff - used by C++ to read geotiff
libtiff - libgeotiff depends on this.
gcc-4.5 or higher, for the -std=c++0x option. Will not compile without it.
py27-yaml - only used by timing.py to save timing results.
google-perftools - optional for doing heap profiling

Installation:

Running "make all" will build three things: a test pure C++ program from main.cpp,
a libraster_stats.so which has the C++ implementation of union-find, and 
raster_stats.so, which is the Python native wrapper created by Boost.Python.

So far, this has been built under OS X 10.7 and Ubuntu.
On Ubuntu
CXX=g++-4.5
Geotiff.h is installed in /opt/include/geotiff/geotiff.h instead
of out in /opt/include, like it is on the mac.

On Mac
CXX=g++-mp-4.5
You have to explicitly link in the Python library because
libboost_python doesn't carry it along.


Testing:

Run unit tests with:

  python traster.py

One of these will fail. It's OK.

Timing:

See timing help with "python timing.py --help". It has two modes, either
time the functions itself or just run so you can do profiling from outside
with time or other performance tools. The latter mode is run with "--heap".

