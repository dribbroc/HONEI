AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

lib_LTLIBRARIES = libhoneibackendscuda.la

NVCC = nvcc
NVCCFLAGS = -g --ptxas-options=-v -O3 -I$(top_srcdir) $(DEBUGDEF) $(PROFILERDEF)

.cu.lo: libtool-hack.in
	$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c -o $(patsubst %.lo,%.o,$@) $<
	sed -e 's/@NAME@/$(patsubst %.lo,%.o,$@)/' < libtool-hack.in > $@

memory_backend_TEST_SOURCES = memory_backend_TEST.cc
memory_backend_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendscuda.la \
	$(DYNAMIC_LD_LIBS)
memory_backend_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

libhoneibackendscuda_la_SOURCES = \
				 cuda_util.hh \
				 operations.hh \
				 difference.cu \
				 dot_product.cu \
				 element_inverse.cu \
				 element_product.cu \
				 product.cu \
				 scale.cu \
				 scaled_sum.cu \
				 sum.cu \
				 transfer.cu \
				 transfer.hh \
				 memory_backend.cc
libhoneibackendscuda_la_LIBADD = \
				 $(top_builddir)/honei/util/libhoneiutil.la \
				-lcudart

libhoneibackendscuda_includedir = $(includedir)/honei/backends/cuda/

libhoneibackendscuda_include_HEADERS = operations.hh cuda_util.hh transfer.hh

TESTS =  memory_backend_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh

check_PROGRAMS = $(TESTS)

