CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF)

noinst_LIBRARIES = libla_spe.a
libla_spe_a_SOURCES = \
	dot_product_dense_dense_float.cc \
	element_product_dense_dense_float.cc \
	product_dense_matrix_dense_matrix_float.cc \
	sum_dense_dense_float.cc \
	sum_dense_scalar_float.cc \
	difference_dense_dense_float.cc \
	product_dense_matrix_dense_vector_float.cc \
	product_banded_matrix_dense_vector_float.cc \
	sum_dense_sparse_float.cc \
	scale_dense_float.cc \
	reduction_sum_dense_float.cc \
	reduction_max_dense_float.cc \
	reduction_min_dense_float.cc \
	norm_max_dense_float.cc \
	norm_l_one_dense_float.cc \
	norm_l_two_dense_float.cc \
	element_inverse_float.cc \
	scaled_sum_dense_dense_float.cc
libla_spe_a_LIBADD = \
	$(top_srcdir)/honei/cell/libutil/libutil_spe.a

libla_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
