AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(OPENCLDEF) \
	$(ITANIUMDEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(CUBLASDEF) \
	$(DEBUGDEF) \
	$(BOOSTDEF) \
	$(PROFILERDEF) \
	$(DUNEDEF) \
	$(MPIDEF) \
	$(GMPDEF) \
	-DHONEI_SOURCEDIR='"$(top_srcdir)"' \
	-DHONEI_BUILDDIR='"$(top_builddir)"'

BACKEND_LIBS = \
       $(top_builddir)/honei/backends/multicore/libhoneibackendsmulticore.la

if CELL

CELLFILES = celllist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/cell/ppe/libhoneibackendscellppe.la \
	$(top_builddir)/honei/backends/cell/spe/libhoneibackendscellspe.la

endif

if CUDA

CUDAFILES = cudalist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/cuda/libhoneibackendscuda.la

endif

if SSE

SSEFILES = sselist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/sse/libhoneibackendssse.la

endif

if ITANIUM

ITANIUMFILES = itaniumlist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/itanium/libhoneibackendsitanium.la

endif

if OPENCL

OPENCLFILES = opencllist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/opencl/libhoneibackendsopencl.la

endif

if MPI

BACKEND_LIBS += \
	$(top_builddir)/honei/backends/mpi/libhoneibackendsmpi.la

endif

OSKILIB = /home/user/dribbroc/nobackup/poski/lib/oski
OSKILIB_SHARED = -I/home/user/dribbroc/nobackup/oski/include -Wl,-rpath -Wl,$(OSKILIB) -L$(OSKILIB) \
		 `cat $(OSKILIB)/site-modules-shared.txt` -loski
POSKILIB = /home/user/dribbroc/nobackup/poski/lib/poski
POSKILIB_SHARED = -I/home/user/dribbroc/nobackup/poski/include -Wl,-rpath -Wl,$(POSKILIB) -L$(POSKILIB) \
		  -lposki
bin_PROGRAMS = poski

poski_SOURCES = poski.cc
poski_CXXFLAGS = -I/home/user/dribbroc/nobackup/poski/include -fopenmp
poski_LDFLAGS = \
	$(OSKILIB_SHARED) $(POSKILIB_SHARED) -lm -fopenmp
poski_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/math/libhoneimath.la \
	$(BACKEND_LIBS)

