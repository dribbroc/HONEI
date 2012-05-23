AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(OPENCLDEF) \
	$(CUDA_DOUBLEDEF) \
	$(CUBLASDEF) \
	$(MPIDEF) \
	$(GMPDEF) \
	$(DEBUGDEF) \
	$(BOOSTDEF) \
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

if MPI

SSEFILES = mpilist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/mpi/libhoneibackendsmpi.la

endif

bin_PROGRAMS = mpi-test honei-mpi honei-mpi-ring honei-mpi-pollutant

mpi_test_SOURCES = mpi_test.cc mpi_test.hh
mpi_test_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)

honei_mpi_SOURCES = honei_mpi.cc honei_mpi.hh
honei_mpi_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)

honei_mpi_ring_SOURCES = honei_mpi_ring.cc honei_mpi_ring.hh
honei_mpi_ring_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)

honei_mpi_pollutant_SOURCES = honei_mpi_pollutant.cc honei_mpi_pollutant.hh
honei_mpi_pollutant_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)
