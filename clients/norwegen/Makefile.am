AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(CUBLASDEF) \
	$(MPIDEF) \
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

if MPI

SSEFILES = mpilist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/mpi/libhoneibackendsmpi.la

endif

bin_PROGRAMS = norwegen


norwegendir = $(top_srcdir)/clients/norwegen
#norwegen_HEADERS = honei_mpi_ring.hh
norwegen_SOURCES = norwegen.cc
norwegen_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(top_builddir)/honei/math/libhoneimath.la \
	$(BACKEND_LIBS)
