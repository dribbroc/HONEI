HONEI README FILE
-----------------


INSTALLING
----------

* Decompress and untar the honei archive.
* Run './configure' in the untarred honei archive's top-level directory.
* Run 'make' to build all libraries and applications.
* Optionally run 'make check TYPE=quick'.

To adapt the honei installation to your system's features, you can pass
several parameters to the configure script.
A full list can be obtained with './configure --help'.

The most common parameters are (see below for dependencies and additional
details):
--enable-cell           Build Cell BE support
--enable-cuda           Build CUDA support
--enable-cuda_double    Add CUDA support with double precission
--enable-sse            Build SSE support
--enable-opencl         Build OpenCL support
--enable-mpi            Build MPI support
--with-cublas           Use the nvidia CuBlas
--with-boost            Use boost instead of tr1 headers
--with-visual           Include GUI support (required for the GUI applications)

Currently only common versions of the GNU Compiler Collection's
c++ compiler are fully supported and tested. Examples include
gcc releases 4.1.x (e.g. OpenSUSE Linux 10.2),
gcc releases 4.2.x (e.g. Ubuntu 8.04),
gcc releases 4.3.x (e.g. Kubuntu 8.10)
Consequently, only compiler flags for GCC are discussed in this file.


COMPILATION DEPENDENCIES AND ADDITIONAL DETAILS
-----------------------------------------------

Compiling individual backends typically requires more fine-tuning at
the configure stage than just passing --enable-XYZ:

- Optimised compiler settings:
  To create an optimised built for your machine (-march=, -mtune=, -ffoobar),
  pass these flags in CXXFLAGS to configure
- GUI and visualisation
  GUIs for the applications honei-poisson and honei-swe need OpenGL, X11 and (free-)glut
  headers and libraries. The applications are only built with --with-visual.
  Necessary files are installed with your graphics card driver, and freeglut is
  part of most Linux distributions, so we assume they have been set up properly
  prior to building HONEI.
- SSE backend:
  You need a machine with at least SSE2 support, such as AMD Athlon
  or Intel Pentium 4 or better.
  You must add the following gcc option to configure:
  CXXFLAGS="-msse2 -mfpmath=sse"
- CUDA backend:
  You need a machine with a CUDA-enabled GPU, and you need to use
  a Linux distribution supported by NVIDIA, see
          http://www.nvidia.com/cuda
  for details. Download and install CUDA and the corresponding device
  driver.
  Set the following environment variables according to your installation,
  ideally in your .shellrc configuration file. We assume here CUDA has been
  installed to /usr/local/cuda/2.0:
          $CUDA_INC_PATH -> /usr/local/cuda/2.0/include/
          $CUDA_LIB_PATH -> /usr/local/cuda/2.0/lib/
          $PATH          -> $PATH:/usr/local/cuda/2.0/bin/
  When configuring, add "-I$CUDA_INC_PATH" to CXXFLAGS,
  and "-L$CUDA_LIB_PATH" to LDFLAGS.
  and 'make' as usual.
- CELL Backend:
  The Cell backend relies heavily on the IBM Cell SDK 2.1.
  Therefore the IBM Cell SDK (at least version 2.1) and libspe2 should be installed properly,
  see http://www-128.ibm.com/developerworks/power/cell/ and
  http://www.bsc.es/plantillaH.php?cat_id=103 for details.
  When configuring, add CC=ppu-gcc CXX=ppu-g++ LD=ppu-ld
  and 'make' as usual.

- Examples
  SSE only, all applications, tuned for Core2Duo: ./configure --enable-sse --with-visual CXXFLAGS="-O3 -msse3 -mfpmath=sse -march=nocona"

  SSE+CUDA, all applications, "generic SSE": ./configure --enable-sse --enable-cuda --with-visual CXXFLAGS="-msse2 -mfpmath=sse -I$CUDA_INC_PATH" LDFLAGS="-L$CUDA_LIB_PATH"

  Cell only: ./configure --enable-cell CXXFLAGS="-O2 -fno-inline" CC=ppu-gcc CXX=ppu-g++ LD=ppu-ld


APPLICATIONS
------------
For building and using the GUI applications (note, that they require '--with-visual'), as well as running tests and benchmarks,
we refer to the associated README files in the following folders:

FOLDER                              DESCRIPTION
benchmark/                          HONEI benchmarks are used to create timing and/or transfer rate and/or performance information.
                                    In particular, this folder contains all benchmarks used throughout the paper.
clients/poisson                     A GUI application to visualise results computed by the HONEI multigrid solver.
clients/swe/                        A GUI application to visualise results computed by the HONEI SWE solver.
honei/math/                         The HONEI math library folder. It also contains solvers and tests.
honei/swe/                          The HONEI swe  library folder. It also contains solvers and tests.


RUNTIME PARAMETERS
------------------

To tune HONEI to your systems specifications, you can configure runtime parameters,
e.g. the number of SPEs used when executing Cell programs.
The easiest way to do so is to copy the file honei/honeirc from the honei tarball
to your home directory under the name .honeirc. Now you can edit ~/.honeirc and
adapt all settings to your needs.
