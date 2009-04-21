CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS= $(SPU_CXXFLAGS) -I$(top_srcdir)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = liblbm_spe.a
liblbm_spe_a_SOURCES = \
	eq_dist_grid_float.cc \
	extraction_grid_float.cc \
	operations.hh

liblbm_spe_a_LIBADD = \
	$(top_builddir)/honei/backends/cell/spe/libutil/libutil_spe.a

liblbm_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti
