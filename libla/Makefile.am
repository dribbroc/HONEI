CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir)

lib_LTLIBRARIES = libla.la

libla_la_SOURCES = \
	banded_matrix.hh \
	dense_matrix.hh \
	dense_vector.hh \
	element_iterator.hh \
	matrix.hh \
	scalar_product.cc scalar_product.hh \
	tags.hh \
	vector.hh \
	vector_error.cc vector_error.hh \
	vector_iterator.hh \
	vector_norm.hh
libla_la_LIBADD = \
	$(top_builddir)/libutil/libutil.la

pg512_includedir = $(includedir)/pg512/
pg512_include_HEADERS = \
	banded_matrix.hh \
	dense_matrix.hh \
	dense_vector.hh \
	element_iterator.hh \
	matrix.hh \
	scalar_product.hh \
	tags.hh \
	vector.hh \
	vector_error.hh \
	vector_iterator.hh \
	vector_norm.hh
