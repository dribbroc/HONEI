AM_CXXFLAGS = -I$(top_srcdir)

MAINTAINERCLEANFILES = Makefile.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
CLEANFILES = *~ \
	     *.outcl

DEFS = \
	$(OPENCLDEF) \
	$(DEBUGDEF) \
	$(BOOSTEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(PROFILERDEF) \
	-DHONEI_SOURCEDIR='"$(abs_top_srcdir)"'

CL_OUT = $(CL_SOURCEN:.cl=.outcl)

all-local: $(CL_OUT)

%.outcl:%.cl oclc
	@$(top_builddir)/honei/backends/opencl/oclc $<
	@echo 'Mnemonic for autotools' > $(<:.cl=.outcl)


lib_LTLIBRARIES = libhoneibackendsopencl.la

libhoneibackendsopencl_la_LIBADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	-lOpenCL

libhoneibackendsopencl_la_SOURCES = defect.cc \
				    difference.cc \
				    dot_product.cc \
				    element_product.cc \
				    fill.cc \
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
CL_SOURCEN =	                    defect.cl \
				    difference.cl \
				    dot_product.cl \
				    element_product.cl \
				    fill.cl \
				    norm.cl \
				    scaled_sum.cl \
				    scale.cl \
				    sum.cl \
				    product.cl




libhoneibackendsopencl_includedir = $(includedir)/honei/backends/opencl/
libhoneibackendsopencl_include_HEADERS =  operations.hh opencl_backend.hh opencl_backend-impl.hh

bin_PROGRAMS = oclc
oclc_SOURCES = oclc.cc
oclc_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	-lOpenCL \
	$(top_builddir)/honei/backends/opencl/libhoneibackendsopencl.la

opencl_backend_TEST_SOURCES = opencl_backend_TEST.cc
opencl_backend_TEST_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendsopencl.la \
	-lOpenCL \
	$(DYNAMIC_LD_LIBS)
opencl_backend__TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

memory_backend_TEST_SOURCES = memory_backend_TEST.cc
memory_backend_TEST_LDADD = \
	$(top_builddir)/honei/util/libhoneiutil.la \
	libhoneibackendsopencl.la \
	-lOpenCL \
	$(DYNAMIC_LD_LIBS)
memory_backend__TEST_CXXFLAGS = -I$(top_srcdir) $(AM_CXXFLAGS)

TESTS =  opencl_backend_TEST memory_backend_TEST
TESTS_ENVIRONMENT = env BACKENDS="$(BACKENDS)" TYPE=$(TYPE) bash $(top_srcdir)/honei/util/run.sh

check_PROGRAMS = $(TESTS)
