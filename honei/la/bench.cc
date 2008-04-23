#define ALLBENCH
#include <benchmark/benchmark.cc> 
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.cc>

#include <string>
#include <tr1/memory>

#include <honei/la/dot_product_BENCHMARK.cc>
// #include <la/vector_scaled_sum_BENCHMARK.cc>
#include <honei/la/element_inverse_BENCHMARK.cc>
#include <honei/la/element_product_BENCHMARK.cc>
#include <honei/la/product_BENCHMARK.cc>
#include <honei/la/reduction_BENCHMARK.cc>
#include <honei/la/scale_BENCHMARK.cc>
#include <honei/la/scaled_sum_BENCHMARK.cc>
#include <honei/la/sum_BENCHMARK.cc>
#include <honei/la/difference_BENCHMARK.cc>
