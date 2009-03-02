CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS= $(SPU_CXXFLAGS) -I$(top_srcdir)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = libmath_spe.a
libmath_spe_a_SOURCES = \
	scaled_product_sum_norm_float.cc \
	sqrt_dense_float.cc \
	operations.hh

libmath_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti
