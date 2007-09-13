#define ALLBENCH
#include <benchmark/benchmark.cc> 
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>

#include <string>
#include <tr1/memory>

#include <libla/dot_product_BENCHMARK.cc>
// #include <libla/vector_scaled_sum_BENCHMARK.cc>
#include <libla/element_inverse_BENCHMARK.cc>
#include <libla/element_product_BENCHMARK.cc>
#include <libla/product_BENCHMARK.cc>
#include <libla/reduction_BENCHMARK.cc>
#include <libla/scale_BENCHMARK.cc>
#include <libla/scaled_sum_BENCHMARK.cc>
#include <libla/sum_BENCHMARK.cc>
