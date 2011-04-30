AM_CXXFLAGS = -I$(top_srcdir)

MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = defect.cl \
	     difference.cl \
	     dot_product.cl \
	     element_product.cl \
	     norm.cl \
	     scaled_sum.cl \
	     scale.cl \
	     sum.cl \
	     product.cl
CLEANFILES = *~

DEFS = \
	$(OPENCLDEF) \
	$(DEBUGDEF) \
	$(BOOSTEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(PROFILERDEF) \
	-DHONEI_SOURCEDIR='"$(abs_top_srcdir)"'

lib_LTLIBRARIES = libhoneibackendsopencl.la

libhoneibackendsopencl_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	-lOpenCL

libhoneibackendsopencl_la_SOURCES = defect.cc \
				    difference.cc \
				    dot_product.cc \
				    element_product.cc \
				    norm.cc \
				    operations.hh \
				    opencl_backend.cc \
				    opencl_backend-impl.hh \
				    memory_backend.cc \
				    opencl_backend.hh \
				    scale.cc \
				    scaled_sum.cc \
				    sum.cc \
				    product.cc


libhoneibackendsopencl_includedir = $(includedir)/honei/backends/opencl/
libhoneibackendsopencl_include_HEADERS =  operations.hh opencl_backend.hh opencl_backend-impl.hh

opencl_backend_TEST_SOURCES = opencl_backend_TEST.cc
opencl_backend_TEST_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendsopencl.la \
	$(DYNAMIC_LD_LIBS)
opencl_backend__TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

memory_backend_TEST_SOURCES = memory_backend_TEST.cc
memory_backend_TEST_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendsopencl.la \
	$(DYNAMIC_LD_LIBS)
memory_backend__TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

TESTS =  opencl_backend_TEST memory_backend_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/honei/util/run.sh

check_PROGRAMS = $(TESTS)
