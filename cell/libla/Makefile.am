CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = spu-g++

noinst_LIBRARIES = libla_spe.a
libla_spe_a_SOURCES = \
	dense_dense_float_dot_product.cc \
	dense_dense_float_element_product.cc \
	dense_dense_float_matrix_product.cc \
	dense_dense_float_sum.cc \
	dense_dense_float_matrix_vector_product.cc \
	dense_sparse_float_sum.cc \
    dense_float_vector_reduction_sum.cc
libla_spe_a_CXXFLAGS = -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
