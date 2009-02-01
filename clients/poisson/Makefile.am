AM_CXXFLAGS = -I$(top_srcdir)

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

bin_PROGRAMS = honei-poisson
honei_poisson_SOURCES = \
			honei_poisson.hh \
			honei_poisson.cc \
			scenario_controller.hh \
			scenario_controller_base.hh
honei_poisson_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/math/libhoneimath.la \
	$(BACKEND_LIBS) \
	-lGL \
	-lglut \
	-lGLU
