#bin_PROGRAMS	= spai convert

#spai_SOURCES		= command_line.c test_spai.c basics.h
#spai_LDADD		= -L. -lspai $(top_builddir)/honei/math/spai/lapack/liblapack.a $(top_builddir)/honei/math/spai/blas/libblas.a -lc -lstdc++ -lgcc_s -ldl -lgfortran
#@FLIBS@ @LIBS@ @C_F77LINK_MAGIC@
#spai_DEPENDENCIES	= libspai.a

#convert_SOURCES		= convert.c
#convert_LDADD		= -L. -lspai @LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@ @LIBS@ @C_F77LINK_MAGIC@
#convert_LDADD		= -L. -lspai $(top_builddir)/honei/math/spai/lapack/liblapack.a $(top_builddir)/honei/math/spai/blas/libblas.a -lc -lstdc++ -lgcc_s -ldl -lgfortran
#convert_DEPENDENCIES	= libspai.a


lib_LTLIBRARIES	= libspai.la

libspai_la_SOURCES	= basics.c bicgstab.c blocks.c com_server.c		\
  debug.c hash_table.c index_set.c						\
  load_balance.c matrix.c mmio.c qr.c read_mm_matrix.c spai.c spai_sub.c	\
  timing.c variable_blocks.c vector.c						\
										\
  basics.h bicgstab.h config.h							\
  blocks.h com_server.h command_line.h debug.h hash_table.h			\
  index_set.h load_balance.h matrix.h mmio.h qr.h read_mm_matrix.h		\
  spai.h spai_error.h spai_sub.h timing.h variable_blocks.h vector.h


# fortran.h


CLEANFILES = test_spai_01.check test_spai_par_01.check

#TESTS      = @SRCTESTS@
#EXTRA_DIST = test_spai_01.sh test_spai_01.out test_spai_par_01.sh.in	\
EXTRA_DIST = local_lapack.h
