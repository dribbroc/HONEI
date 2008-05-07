CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
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
	$(top_builddir)/honei/backends/cuda/libhoneibackendscuda.so

endif

if SSE

SSEFILES = sselist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/sse/libhoneibackendssse.la

endif

bin_PROGRAMS = honei-swe
honei_swedir = $(top_srcdir)/clients/swe
honei_swe_HEADERS = scenario_controller.hh
honei_swe_SOURCES = honei_swe.cc
honei_swe_LDADD = \
	$(top_srcdir)/honei/util/libhoneiutil.la \
	$(top_srcdir)/honei/la/libhoneila.la \
	$(top_srcdir)/honei/swe/libhoneiswe.la \
	$(BACKEND_LIBS) \
	-lGL \
	-lglut \
	-lGLU
