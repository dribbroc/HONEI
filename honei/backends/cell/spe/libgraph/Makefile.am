CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = libgraph_spe.a
libgraph_spe_a_SOURCES = \
	node_distance_float_wfr.cc \
	node_distance_float_wkk.cc
libgraph_spe_a_LIBADD = \
	$(top_srcdir)/honei/backends/cell/spe/libutil/libutil_spe.a

libgraph_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
