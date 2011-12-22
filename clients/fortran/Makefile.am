AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~ feat
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(DEBUGDEF) \
	$(BOOSTDEF) \
	$(GMPDEF) \
	$(PROFILERDEF)

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

if OPENCL

OPENCLFILES = opencllist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/opencl/libhoneibackendsopencl.la

endif

lib_LTLIBRARIES = libhoneifortran.la

libhoneifortran_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/math/libhoneimath.la \
	$(BACKEND_LIBS)

libhoneifortran_la_SOURCES = ftoh.cc

#libhoneifortran_includedir = $(includedir)/honei/backends/sse/
#libhoneifortran_include_HEADERS =

feat: libhoneifortran.la
	gfortran $(CXXFLAGS) -o feat feat.f90 -lhoneifortran -L.libs
	@echo To run the test application, run
	@echo LD_LIBRARY_PATH=\"{DOLLAR}LD_LIBRARY_PATH:.libs\" ./feat
