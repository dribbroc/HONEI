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
	$(top_builddir)/honei/backends/cuda/libhoneibackendscuda.la

endif

if SSE

SSEFILES = sselist
BACKEND_LIBS += \
	$(top_builddir)/honei/backends/sse/libhoneibackendssse.la

endif

bin_PROGRAMS = honei-lbm
honei_lbmdir = $(top_srcdir)/clients/lbm
honei_lbm_HEADERS = scenario_controller.hh
honei_lbm_SOURCES = honei_lbm.cc
honei_lbm_LDADD = \
	$(top_srcdir)/honei/util/libhoneiutil.la \
	$(top_srcdir)/honei/la/libhoneila.la \
	$(top_srcdir)/honei/swe/libhoneiswe.la \
	$(BACKEND_LIBS) \
	-lGL \
	-lglut \
	-lGLU
