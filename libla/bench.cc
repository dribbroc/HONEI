#define ALLBENCH
#include <benchmark/benchmark.cc> 
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>

#include <string>
#include <tr1/memory>

#include <libla/scalar_product_BENCHMARK.cc>
#include <libla/vector_scaled_sum_BENCHMARK.cc>
#include <libla/matrix_element_inverse_BENCHMARK.cc>
#include <libla/matrix_elementwise_product_BENCHMARK.cc>
#include <libla/matrix_product_BENCHMARK.cc>
#include <libla/matrix_row_sum_vector_BENCHMARK.cc>
#include <libla/scalar_matrix_product_BENCHMARK.cc>
#include <libla/scalar_matrix_sum_BENCHMARK.cc>
