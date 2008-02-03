#define ALLBENCH
#include <benchmark/benchmark.cc> 
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/matrix_error.cc>

#include <string>
#include <tr1/memory>

#include <honei/libla/dot_product_BENCHMARK.cc>
// #include <libla/vector_scaled_sum_BENCHMARK.cc>
#include <honei/libla/element_inverse_BENCHMARK.cc>
#include <honei/libla/element_product_BENCHMARK.cc>
#include <honei/libla/product_BENCHMARK.cc>
#include <honei/libla/reduction_BENCHMARK.cc>
#include <honei/libla/scale_BENCHMARK.cc>
#include <honei/libla/scaled_sum_BENCHMARK.cc>
#include <honei/libla/sum_BENCHMARK.cc>
#include <honei/libla/difference_BENCHMARK.cc>
