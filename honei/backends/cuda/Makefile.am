AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2

lib_LTLIBRARIES = libhoneibackendscuda.la

NVCC = nvcc
NVCCFLAGS = -g

.cu.lo: libtool-hack.in
	$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c -o $(patsubst %.lo,%.o,$@) $<
	sed -e 's/@NAME@/$(patsubst %.lo,%.o,$@)/' < libtool-hack.in > $@

libhoneibackendscuda_la_SOURCES = \
				 operations.hh \
				 element_inverse.cu \
				 element_product.cu \
				 scale.cu \
				 scaled_sum.cu \
				 sum.cu
libhoneibackendscuda_la_LIBADD = \
				-lcudart

libhoneibackendscuda_includedir = $(includedir)/honei/backends/cuda/

libhoneibackendscuda_include_HEADERS = operations.hh
