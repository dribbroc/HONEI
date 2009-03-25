AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(ITANIUMDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

lib_LTLIBRARIES = libhoneibackendsitanium.la

libhoneibackendsitanium_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendsitanium_la_SOURCES = operations.hh \
				 scaled_sum.cc

libhoneibackendsitanium_includedir = $(includedir)/honei/backends/itanium/
libhoneibackendsitanium_include_HEADERS =  operations.hh
