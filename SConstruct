'''
This builds the raster_measure library. It's an SCons configuration file.

Build: "scons"
Targets:
  - cpp, builds all the pure C++ (Default)
  - wrapper, builds the Python wrapper dll
  - tiff_test, cluster_test are two test cases.
Clean: "scons -c"
Reset cache: "rm -rf .scon*"
Debug: "scons --echo"

There is a configuration file called default.cfg. If you want to change
values for a particular OS, like "linux2" or "darwin", then edit
linux2.cfg or darwin.cfg. If you have defaults just for your machine,
then create local.cfg and put those defaults there.

Want to know which of the variables in the cfg files were read? Look
at used.cfg after the build is done.
'''

import os
import sys
import subprocess
import datetime
import sysconfig
import SConsDebug
import SConsCheck
from SConsVars import cfg
import logging


logger=logging.getLogger('SconsRoot')
logging.basicConfig(level=logging.DEBUG)

# The goals of the build decide what we find and include.
# 1. Command line determines what targets we need.
# 2. Targets have requirements.
# 3. Availability can determine optional things to include (like tcmalloc).

# Set preferred directories here. There is a hierarchy of what gets chosen.
# 1. Set from command line
# 2. Set from environment variable.
# 3. Set for this particular machine in some config file.
# 4. Set for this architecture.
# 5. Set as common default, such as /usr/include.


# Command-line Options that are arguments to SCons.
file_cpu=int(cfg.get('General','num_jobs'))
num_cpu=int(os.environ.get('NUM_CPU',file_cpu))
SetOption('num_jobs', num_cpu)
logger.info('running with -j %d' % int(GetOption('num_jobs')))

AddOption('--echo', dest='echo', action='store_true',default=False,
          help='This echoes what scons runs as tests in order to debug failed builds.')
AddOption('--tcmalloc', dest='tcmalloc', type='int',action='store',default=0,
      help='Whether to use tcmalloc to debug the heap.')
AddOption('--cpp_compiler', dest='cpp_compiler', action='store', default=None,
          help='Set the C++ compiler.')
AddOption('--c_compiler', dest='c_compiler', action='store', default=None,
          help='Set the C compiler.')

env=Environment()

if GetOption('echo'):
    env['SPAWN']=SConsDebug.echospawn
cpp_compiler=GetOption('cpp_compiler')
cpp_compiler=cpp_compiler or cfg.get('General','cpp_compiler')
cpp_compiler=cpp_compiler or SConsCheck.latest_gcc(env)

c_compiler=GetOption('c_compiler')
c_compiler=c_compiler or cfg.get('General','c_compiler')

if cpp_compiler:
    env['CXX']=cpp_compiler
if c_compiler:
    env['CC']=cpp_compiler


def DisallowSubst():
    #AllowSubstExceptions()
    # Cannot use AllowSubstExceptions without defining these variables
    # which the default tool leaves undefined.
    env.SetDefault(CPPFLAGS=[])
    env.SetDefault(CPPDEFINES=[])
    env.SetDefault(CXXCOMSTR='')
    env.SetDefault(RPATH=[])
    env.SetDefault(LIBPATH=[])
    env.SetDefault(LIBS=[])
    env.SetDefault(LINKCOMSTR='')
    env.SetDefault(SHCXXCOMSTR='')
    env.SetDefault(SHLINKCOMSTR='')



DisallowSubst()

env.AppendUnique(CCFLAGS=['-fPIC'])
env.AppendUnique(LINKFLAGS=['-fPIC'])
env.AppendUnique(CPPPATH=['.'])
optimization=cfg.get('General','optimization').strip().split()
if optimization:
    env.AppendUnique(CCFLAGS   = optimization )
    env.AppendUnique(LINKFLAGS = optimization )


tiffmaker=SConsCheck.which('pnmtotiff')
if tiffmaker:
    logger.debug('Creating a builder for TIFF from PGM.')
    MakeTIFF = Builder(action = '%s $SOURCE > $TARGET' % tiffmaker)
else:
    logger.warn('The tiff_test input file needs you to install netpbm '+
                'for pnmtotiff. You can run other tests without this.')
    MakeTIFF=None

# These are the boost libraries needed by the code.
boost_libs  = ['boost_system-mt','boost_chrono-mt','boost_random-mt',
    'boost_program_options-mt','boost_filesystem-mt']


conf = Configure(env, custom_tests = {
     'CheckBoost' : SConsCheck.CheckBoost(boost_libs),
     'CheckCPP11' : SConsCheck.CheckCPP11(),
     'CheckGeoTIFF' : SConsCheck.CheckGeoTIFF(),
     'CheckTBB' : SConsCheck.CheckTBB()
    })

if GetOption('echo'):
    conf.logstream = sys.stdout

if not conf.CheckCXX():
    logger.debug('CXX %s' % env['CXX'])
    logger.debug('CXXCOM %s' % env['CXXCOM'])
    logger.error('The compiler isn\'t compiling.')
    Exit(1)

failure_cnt=0

if not conf.CheckCPP11():
    logger.error('Compiler does not support modern C++ standards (C++23/C++20/C++17).')
    failure_cnt+=1

if not conf.CheckBoost():
    logger.error('Boost not found')
    failure_cnt+=1

if not conf.CheckGeoTIFF():
    logger.error('GeoTIFF not found')
    failure_cnt+=1

if not conf.CheckLib('m',language='C++'):
    logger.error('No math library?!')
    failure_cnt+=1

if not conf.CheckLib('pthread',language='C++'):
    logger.error('We need pthread.')
    failure_cnt+=1

if not conf.CheckLib('libhdf5',language='C'):
    logger.error('Cannot load HDF5 library')
    failure_cnt+=1

if not conf.CheckLib('libhdf5_cpp',language='C'):
    logger.error('Cannot load HDF5 library')
    failure_cnt+=1

if not conf.CheckLib('libhdf5_hl',language='C'):
    logger.error('Cannot load HDF5 library')
    failure_cnt+=1

if not conf.CheckLib('libhdf5_hl_cpp',language='C'):
    logger.error('Cannot load HDF5 library')
    failure_cnt+=1

if conf.CheckTBB():
    # Not necessarily a failure.
    tbb_exists=True
else:
    tbb_exists=False
    
if failure_cnt:
    Exit(2)

env = conf.Finish()


# The Python environment is for building the Python wrappers.
py_env = env.Clone()

py_conf = Configure(py_env, custom_tests = {
     'CheckPython' : SConsCheck.CheckPython(),
     'CheckBoostPython' : SConsCheck.CheckBoostPython(),
     'CheckNumpy' : SConsCheck.CheckNumpy()
     })

if GetOption('echo'):
    py_conf.logstream = sys.stdout

failure_cnt = 0

if not py_conf.CheckPython():
    logger.error('Could not find Python development directories.')
    failure_cnt+=1

if not py_conf.CheckBoostPython():
    logger.error('Could not find Boost Python development directories.')
    failure_cnt+=1

if not py_conf.CheckNumpy():
    logger.error('Could not find Numpy development directories.')
    failure_cnt+=1

if failure_cnt:
    cfg.write('used.cfg')
    #Exit(3)
    
py_env = py_conf.Finish()

# Testing build environment
test_env = env.Clone()
# Use dynamic library variant of the unit testing framework.
test_env.Append(CCFLAGS = ['-DBOOST_TEST_DYN_LINK'])

if tiffmaker:
    test_env.Append(BUILDERS = {'Tif' : MakeTIFF})

test_conf = Configure(test_env, custom_tests = {})
if GetOption('echo'):
    test_conf.logstream = sys.stdout
test_conf.CheckLib('boost_unit_test_framework',language='C++')
test_env = test_conf.Finish()


# Now begin building.
common = ['io_geotiff.cpp','cluster.cpp','io_ppm.cpp','timing.cpp',
          'timing_harness.cpp', 'cluster_generic.cpp']
if tbb_exists:
    common += ['cluster_tbb.cpp']

stats = env.SharedLibrary(target='raster_stats',
                          source=common)
                          
python_wrap_sources=['raster_wrap.cpp']
wrapper=py_env.SharedLibrary(target='raster_wrap', source=[stats]+python_wrap_sources)
Alias('wrapper',wrapper)

raster = env.Program(target='raster',source=['main.cpp',stats])

cluster_test = test_env.Program('cluster_test',source=['cluster_test.cpp',stats])

if tbb_exists:
    cluster_tbb_test = test_env.Program('cluster_tbb_test',source=['cluster_tbb_test.cpp',stats])

logger.debug('About to add MakeTIFF')
if tiffmaker:
    logger.debug('Adding tif production to the build.')
    feep_tif = test_env.Tif('feep.tif','feep.pgm')
tiff_test = test_env.Program('tiff_test',source=['io_geotiff_test.cpp',stats])
mac_main = test_env.Program('mac_main',source=['mac_main.cpp',stats])

mac_test = test_env.Program('mac_test',source=['mac_test.cpp'])

# This alias covers several environments, so maybe it should be attached to
# the implicit environment.
cpp_includes = [tiff_test, cluster_test, raster,mac_main]
if tiffmaker:
    cpp_includes.append(feep_tif)
if tbb_exists:
    cpp_includes.append(cluster_tbb_test)
cpp_target=Alias('cpp', cpp_includes)
Default(cpp_target)
all=Alias('all',cpp_target)

cfg.write('used.cfg')

try:
    svnversion = subprocess.check_output(['svn','info','--xml'])
    if isinstance(svnversion, bytes):
        svnversion = svnversion.decode('utf-8')
except (subprocess.CalledProcessError, FileNotFoundError):
    svnversion = 'No SVN info available'
text_cfg=open('used.cfg','r').read().replace(os.linesep,'\t')
f=open('raster_version.hpp','w')
f.write('''
#define RASTER_STATS_VERSION R"(%s)"

#define RASTER_STATS_CFG R"(%s)"

#define RASTER_STATS_COMPILE_TIME R"(%s)"

''' % (svnversion.replace(os.linesep,''),text_cfg,
       datetime.datetime.now().isoformat()))
f.close()

logger.info('End of first pass of SConstruct')
