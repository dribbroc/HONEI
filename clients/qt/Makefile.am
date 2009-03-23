include $(top_srcdir)/misc/autotroll.mk

AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~ *.moc.cc
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(DEBUGDEF) \
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

bin_PROGRAMS = honei-qt
honei_qt_SOURCES = honei_qt.cc window.hh window.cc gl_widget.hh gl_widget.cc
honei_qt_CXXFLAGS = $(QT_CXXFLAGS) $(QT_INCPATH) $(AM_CXXFLAGS)
honei_qt_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
honei_qt_LDFLAGS  = $(QT_LDFLAGS) $(LDFLAGS)
honei_qt_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	$(top_builddir)/honei/la/libhoneila.la \
	$(top_builddir)/honei/swe/libhoneiswe.la \
	$(top_builddir)/honei/lbm/libhoneilbm.la \
	$(BACKEND_LIBS) \
	$(QT_LIBS) $(LDADD)
nodist_honei_qt_SOURCES = window.moc.cc gl_widget.moc.cc
