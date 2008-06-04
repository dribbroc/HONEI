AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

lib_LTLIBRARIES = libhoneibackendscuda.la

NVCC = nvcc
NVCCFLAGS = -v -g --ptxas-options=-v -O3 -I$(top_srcdir) $(DEBUGDEF) $(PROFILERDEF)

.cu.lo: libtool-hack.in
	$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c -o $(patsubst %.lo,%.o,$@) $<
	sed -e 's/@NAME@/$(patsubst %.lo,%.o,$@)/' < libtool-hack.in > $@

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
				 sum.cu
libhoneibackendscuda_la_LIBADD = \
				 $(top_builddir)/honei/util/libhoneiutil.la \
				-lcudart

libhoneibackendscuda_includedir = $(includedir)/honei/backends/cuda/

libhoneibackendscuda_include_HEADERS = operations.hh cuda_util.hh
