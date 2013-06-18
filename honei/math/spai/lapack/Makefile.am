
lib_LTLIBRARIES = liblapack.la
liblapack_la_SOURCES = lapack_subs.f
if BLAPACK
else
liblapack_la_LIBADD = \
	$(top_builddir)/honei/math/spai/blas/libblas.la
endif
