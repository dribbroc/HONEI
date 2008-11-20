CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS= -I$(top_srcdir) #$(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF) $(PROFILERDEF)

noinst_LIBRARIES = libswe_spe.a
libswe_spe_a_SOURCES = \
	flow_processing_float.cc \
	operations.hh \
	triple.hh

libswe_spe_a_LIBADD = \
	$(top_builddir)/honei/backends/cell/spe/libutil/libutil_spe.a

libswe_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti
