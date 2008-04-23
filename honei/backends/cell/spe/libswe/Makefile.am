CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = #$(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF)

noinst_LIBRARIES = libswe_spe.a
libswe_spe_a_SOURCES = \
	flow_processing_float.cc

libswe_spe_a_LIBADD = \
	$(top_srcdir)/honei/backends/cell/spe/libutil/libutil_spe.a

libswe_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
