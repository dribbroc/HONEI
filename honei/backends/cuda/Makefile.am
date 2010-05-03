AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~ *.linkinfo
MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = libtool-hack.in

lib_LTLIBRARIES = libhoneibackendscuda.la

DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(CUBLASDEF) \
	$(DEBUGDEF) \
	$(PROFILERDEF)

NVCC = nvcc
NVCCFLAGS = -g --ptxas-options=-v -O3 -I$(top_srcdir) $(DEFS)
if CUDA_DOUBLE
  NVCCFLAGS += -arch=sm_13
endif

.cu.lo: libtool-hack.in
	$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c -o $(patsubst %.lo,%.o,$@) $<
	sed -e 's/@NAME@/$(patsubst %.lo,%.o,$@)/' < $(srcdir)/libtool-hack.in > $@

memory_backend_TEST_SOURCES = memory_backend_TEST.cc
memory_backend_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendscuda.la \
	$(DYNAMIC_LD_LIBS)
memory_backend__TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

multi_gpu_TEST_SOURCES = multi_gpu_TEST.cc
multi_gpu_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendscuda.la \
	$(DYNAMIC_LD_LIBS)
multi_gpu_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

gpu_pool_TEST_SOURCES = gpu_pool_TEST.cc
gpu_pool_TEST_LDADD = \
	$(top_builddir)/unittest/libunittest.a \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendscuda.la \
	$(DYNAMIC_LD_LIBS)
gpu_pool_TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)


libhoneibackendscuda_la_SOURCES = \
				 cuda_util.hh \
				 boundary_init_fsi.cu \
				 collide_stream_grid.cu \
				 collide_stream_fsi.cu \
				 defect.cu \
				 difference.cu \
				 dot_product.cu \
				 element_inverse.cu \
				 element_product.cu \
				 extraction_grid.cu \
				 eq_dist_grid.cu \
				 force_grid.cu \
				 memory_backend.cc \
				 multi_gpu.cu \
				 multi_gpu.hh \
				 norm.cu \
				 operations.hh \
				 operations.cu \
				 prolongation.cu \
				 product.cu \
				 restriction.cu \
				 scale.cu \
				 scaled_sum.cu \
				 sum.cu \
				 transfer.cu \
				 transfer.hh \
				 up_vel_dir_grid.cu \
				 gpu_function.cc \
				 gpu_function.hh \
				 gpu_pool.cc \
				 gpu_pool.hh \
				 ticket.cc \
				 ticket.hh
libhoneibackendscuda_la_LIBADD = \
				 $(top_builddir)/honei/util/libhoneiutil.la \
				-lcudart \
				-lcublas

libhoneibackendscuda_includedir = $(includedir)/honei/backends/cuda/

libhoneibackendscuda_include_HEADERS = operations.hh \
				       cuda_util.hh \
				       transfer.hh \
				       multi_gpu.hh \
				       gpu_function.hh \
				       gpu_pool.hh \
				       ticket.hh

TESTS =  memory_backend_TEST multi_gpu_TEST gpu_pool_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/unittest/run.sh

check_PROGRAMS = $(TESTS)

