CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

CXX = $(SPU_CXX)
CXXFLAGS = $(SPU_CXXFLAGS)
DEFS = $(DEBUGDEF)

noinst_LIBRARIES = libla_spe.a
libla_spe_a_SOURCES = \
	dense_dense_float_dot_product.cc \
	dense_dense_float_element_product.cc \
	dense_dense_float_matrix_product.cc \
	dense_dense_float_sum.cc \
	dense_scalar_float_sum.cc \
	dense_dense_float_difference.cc \
	dense_dense_float_matrix_vector_product.cc \
	banded_dense_float_matrix_vector_product.cc \
	dense_sparse_float_sum.cc \
	dense_float_scale.cc \
	dense_float_reduction_sum.cc \
	dense_float_reduction_max.cc \
	dense_float_reduction_min.cc \
	dense_float_norm_max.cc \
	dense_float_norm_l_one.cc \
	dense_float_norm_l_two.cc \
	float_element_inverse.cc \
	dense_dense_float_scaled_sum.cc

libla_spe_a_CXXFLAGS = -Wall -msafe-dma -fno-exceptions -fno-rtti

bench:
quickcheck:
quickcheck-sse:
quickcheck-cell:
quickcheck-mc:
check-sse:
check-cell:
check-mc:
