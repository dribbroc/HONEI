AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

lib_LTLIBRARIES = libhoneibackendssse.la

libhoneibackendssse_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la

libhoneibackendssse_la_SOURCES = operations.hh \
				 difference.cc \
				 dot_product.cc \
				 element_inverse.cc \
				 element_product.cc \
				 norm.cc \
				 reduction.cc \
				 scaled_sum.cc \
				 scale.cc \
				 sum.cc

libhoneibackendssse_includedir = $(includedir)/honei/backends/sse/
libhoneibackendssse_include_HEADERS =  operations.hh
