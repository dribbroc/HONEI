CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

bin_PROGRAMS = honei-swe
honei_swedir = $(top_srcdir)/clients/swe
honei_swe_HEADERS = scenario_controller.hh
honei_swe_SOURCES = honei_swe.cc
honei_swe_LDADD = \
	$(top_srcdir)/honei/util/libhoneiutil.la \
	$(top_srcdir)/honei/libla/libla.la \
	$(top_srcdir)/honei/swe/libhoneiswe.la \
	-lGL \
	-lglut \
	-lGLU
