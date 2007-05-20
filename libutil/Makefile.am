CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir)

lib_LTLIBRARIES = libutil.la

libutil_la_SOURCES = \
	exception.cc exception.hh \
	log.cc log.hh \
	shared_array.cc shared_array.hh \
	stringify.hh

pg512_libutil_includedir = $(includedir)/pg512/libutil
pg512_libutil_include_HEADERS = \
	exception.hh \
	log.hh \
	shared_array.hh \
	stringify.hh
