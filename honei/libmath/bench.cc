#define ALLBENCH
#include <benchmark/benchmark.cc> 
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/matrix_error.cc>

#include <string>
#include <tr1/memory>

#include<poisson_cg_double_dense_BENCHMARK.cc>
#include<poisson_cg_float_banded_BENCHMARK.cc>
#include<poisson_iteref_cg_double_dense_BENCHMARK.cc>
#include<poisson_iteref_cg_float_banded_BENCHMARK.cc>
#include<poisson_iteref_pcg_double_dense_BENCHMARK.cc>
#include<poisson_iteref_pcg_float_banded_BENCHMARK.cc>
#include<poisson_jac_double_dense_BENCHMARK.cc>
#include<poisson_jac_float_banded_BENCHMARK.cc>
#include<poisson_jackernel_double_banded_BENCHMARK.cc>
#include<poisson_jackernel_double_banded_sse_BENCHMARK.cc>
#include<poisson_jackernel_float_banded_BENCHMARK.cc>
#include<poisson_jackernel_float_banded_sse_BENCHMARK.cc>
#include<poisson_pcg_double_dense_BENCHMARK.cc>
#include<poisson_pcg_float_banded_BENCHMARK.cc>
