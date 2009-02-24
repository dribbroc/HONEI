AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = libtool-hack.in

lib_LTLIBRARIES = libhoneibackendscuda.la

NVCC = nvcc
NVCCFLAGS = -g --ptxas-options=-v -O3 -I$(top_srcdir) $(DEBUGDEF) $(PROFILERDEF)

.cu.lo: libtool-hack.in
	$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c -o $(patsubst %.lo,%.o,$@) $<
	sed -e 's/@NAME@/$(patsubst %.lo,%.o,$@)/' < $(srcdir)/libtool-hack.in > $@

memory_backend_TEST_SOURCES = memory_backend_TEST.cc
memory_backend_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendscuda.la \
	$(DYNAMIC_LD_LIBS)
memory_backend_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

libhoneibackendscuda_la_SOURCES = \
				 cuda_util.hh \
				 collide_stream_grid.cu \
				 defect.cu \
				 difference.cu \
				 dot_product.cu \
				 element_inverse.cu \
				 element_product.cu \
				 extraction_grid.cu \
				 eq_dist_grid.cu \
				 force_grid.cu \
				 memory_backend.cc \
				 norm.cu \
				 operations.hh \
				 operations.cu \
				 prolongation.cu \
				 product.cu \
				 restriction.cu \
				 scale.cu \
				 scaled_sum.cu \
				 scaled_product_sum_norm.cu \
				 sum.cu \
				 transfer.cu \
				 transfer.hh \
				 up_vel_dir_grid.cu
libhoneibackendscuda_la_LIBADD = \
				 $(top_builddir)/honei/util/libhoneiutil.la \
				-lcudart

libhoneibackendscuda_includedir = $(includedir)/honei/backends/cuda/

libhoneibackendscuda_include_HEADERS = operations.hh cuda_util.hh transfer.hh

TESTS =  memory_backend_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh

check_PROGRAMS = $(TESTS)

