CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS= $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF)

noinst_LIBRARIES = libutil_spe.a

libutil_spe_a_SOURCES = \
	allocator.cc \
	operations.cc \
	operation_one_float_dma.cc \
	operation_one_float_mail.cc \
	operation_two_float_mail.cc \
	operation_two_float_dma.cc \
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
