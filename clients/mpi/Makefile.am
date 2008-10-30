CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(MPIDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

BACKEND_LIBS =

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

bin_PROGRAMS = mpi-test honei-mpi
mpi_testdir = $(top_srcdir)/clients/mpi
mpi_test_HEADERS = mpi_test.hh
mpi_test_SOURCES = mpi_test.cc
mpi_test_LDADD = \
	$(top_srcdir)/honei/util/libhoneiutil.la \
	$(top_srcdir)/honei/la/libhoneila.la \
	$(top_srcdir)/honei/swe/libhoneiswe.la \
	$(top_srcdir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)

honei_mpidir = $(top_srcdir)/clients/mpi
honei_mpi_HEADERS = honei_mpi.hh
honei_mpi_SOURCES = honei_mpi.cc
honei_mpi_LDADD = \
	$(top_srcdir)/honei/util/libhoneiutil.la \
	$(top_srcdir)/honei/la/libhoneila.la \
	$(top_srcdir)/honei/swe/libhoneiswe.la \
	$(top_srcdir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS)
