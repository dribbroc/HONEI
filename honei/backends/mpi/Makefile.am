AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(MPIDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

lib_LTLIBRARIES = libhoneibackendsmpi.la

libhoneibackendsmpi_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendsmpi_la_SOURCES = operations.hh \
				 operations.cc

libhoneibackendsmpi_includedir = $(includedir)/honei/backends/mpi/
libhoneibackendsmpi_include_HEADERS =  operations.hh
