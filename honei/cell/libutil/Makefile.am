CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS= $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = libutil_spe.a

libutil_spe_a_HEADER = \
	allocator.hh \
	operation_framework.hh \
	profiler.hh \
	transfer.hh
libutil_spe_a_SOURCES = \
	allocator.cc \
	operation_framework.cc \
	operation_one_float_dma.cc \
	operation_one_double_dma.cc \
	operation_one_float_mail.cc \
	operation_one_double_mail.cc \
	operation_two_float_mail.cc \
	operation_two_double_mail.cc \
	operation_two_float_dma.cc \
	operation_two_double_dma.cc \
	operation_three_float_dma.cc \
	profiler.cc \
	transfer.cc

libutil_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
