CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF)

noinst_LIBRARIES = libgraph_spe.a
libgraph_spe_a_SOURCES = \
libgraph_spe_a_LIBADD = \
	$(top_srcdir)/cell/libutil/libutil_spe.a

libgraph_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
