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
				 collide_stream_grid.cc \
				 difference.cc \
				 dot_product.cc \
				 element_inverse.cc \
				 element_product.cc \
				 eq_dist_grid.cc \
				 extraction_grid.cc \
				 norm.cc \
				 product.cc \
				 reduction.cc \
				 scaled_sum.cc \
				 scale.cc \
				 sum.cc \
				 up_vel_dir_grid.cc

libhoneibackendssse_includedir = $(includedir)/honei/backends/sse/
libhoneibackendssse_include_HEADERS =  operations.hh
