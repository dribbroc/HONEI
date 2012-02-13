AM_CXXFLAGS = -I$(top_srcdir)

CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in
DEFS = \
	$(CELLDEF) \
	$(SSEDEF) \
	$(OPENCLDEF) \
	$(ITANIUMDEF) \
	$(CUDADEF) \
	$(CUDA_DOUBLEDEF) \
	$(CUBLASDEF) \
	$(DEBUGDEF) \
	$(BOOSTDEF) \
	$(PROFILERDEF) \
	$(DUNEDEF) \
	$(MPIDEF) \
	$(GMPDEF) \
	-DHONEI_SOURCEDIR='"$(top_srcdir)"' \
	-DHONEI_BUILDDIR='"$(top_builddir)"'

lib_LTLIBRARIES = libbenzi.la

libbenzi_la_SOURCES = sainvsr_driver.F ainvini1.h cd2cs3_f90.h mvucsr1.h rank1n.h reallocr.h rifsr.h srtcs1_f90.h ainvsr2.h \
		      cd2cs_f90.h rank1n2.h realloci.h rifsl3.f.h sparslab.i tvucsr1.h

