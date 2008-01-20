CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = libmath_spe.a
libmath_spe_a_SOURCES = \
	sqrt_dense_float.cc

libmath_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
