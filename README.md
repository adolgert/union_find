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

python3.13 - Modern Python 3 version (Python 3.8+ should work)
numpy and scipy (Python 3 versions)
gdal - used to read geotiff in Python scripts
libgeotiff - used by C++ to read geotiff
libtiff - libgeotiff depends on this
Modern C++ compiler supporting C++23, C++20, or C++17 (GCC 7+, Clang 5+)
  - The build system will automatically detect the best available C++ standard
Boost libraries (with Python 3 support):
  - boost_system, boost_chrono, boost_random, boost_program_options, boost_filesystem
  - boost_python3 (or boost_python313/boost_python-py3)
  - boost_unit_test_framework (for testing)
pyyaml - only used by timing.py to save timing results
google-perftools - optional for doing heap profiling
HDF5 library with C++ bindings
Intel TBB (Threading Building Blocks) - optional for parallel algorithms

Installation:

Running "make all" will build three things: a test pure C++ program from main.cpp,
a libraster_stats.so which has the C++ implementation of union-find, and 
raster_stats.so, which is the Python native wrapper created by Boost.Python.

Tested on modern Linux and macOS systems.

The build system auto-detects:
- The best available C++ compiler (tries modern GCC/Clang first)
- Python 3.13 installation and headers
- Boost.Python library names for Python 3
- Library paths for dependencies

Note: Geotiff.h may be installed in different locations:
- /opt/include/geotiff/geotiff.h
- /usr/include/geotiff.h
- /usr/local/include/geotiff.h
The build system searches common locations automatically.


Testing:

Run unit tests with:

  python3 traster.py

One of these will fail. It's OK.

Timing:

See timing help with "python3 timing.py --help". It has two modes, either
time the functions itself or just run so you can do profiling from outside
with time or other performance tools. The latter mode is run with "--heap".

