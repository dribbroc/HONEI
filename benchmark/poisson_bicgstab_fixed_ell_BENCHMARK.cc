/* vim: set sw=4 sts=4 et foldmethod=syntax : */

//#define SOLVER_VERBOSE 1
#include <honei/math/bi_conjugate_gradients_stabilised.hh>
#include <honei/la/element_inverse.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/transposition.hh>
#include <benchmark/benchmark.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/math/endian_swap.hh>
#include <honei/math/spai.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace std;

template<typename Tag_, typename DT1_, typename Method_>
class BiCGStabPoissonAdvancedEllBENCH
{
};

template <typename Tag_, typename DT1_>
class BiCGStabPoissonAdvancedEllBENCH<Tag_, DT1_, methods::JAC>:
    public Benchmark
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;

    public:
        BiCGStabPoissonAdvancedEllBENCH(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
        }

        virtual void run()
        {
            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";

            std::string A_file(file_base);
            A_file += "A_";
            A_file += stringify(_size);
            A_file += ".ell";
            SparseMatrixELL<DT1_> system(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

            std::string rhs_file(file_base);
            rhs_file += "rhs_";
            rhs_file += stringify(_size);
            DenseVector<DT1_> rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0));

            SparseMatrix<DT1_> tsmatrix3(system.rows(), system.columns());
            for(unsigned long i(0) ; i < system.rows() ; ++i)
            {
                tsmatrix3(i , i) = 0.7 * 1/system(i, i);
            }
            SparseMatrixELL<DT1_> smatrix3(tsmatrix3);

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            unsigned long used_iters(0);

            BENCHMARK(
                      (PBiCGStab<Tag_, methods::VAR>::value(system, rhs, result, smatrix3, 10000ul, used_iters));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q1("PARENG JAC BiCGStab double sse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q1("PARENG JAC BiCGStab double sse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q1("PARENG JAC BiCGStab double sse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_0_q1("PARENG JAC BiCGStab double sse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_0_q2("PARENG JAC BiCGStab double sse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q2("PARENG JAC BiCGStab double sse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q2("PARENG JAC BiCGStab double sse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q2("PARENG JAC BiCGStab double sse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");


BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q1("PARENG JAC BiCGStab double mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q1("PARENG JAC BiCGStab double mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_0_q1("PARENG JAC BiCGStab double mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_0_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q2("PARENG JAC BiCGStab double mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q2("PARENG JAC BiCGStab double mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q2("PARENG JAC BiCGStab double mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q1("PARENG JAC BiCGStab double cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q1("PARENG JAC BiCGStab double cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_0_q1("PARENG JAC BiCGStab double cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_0_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_0_q2("PARENG JAC BiCGStab double cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_0_q2("PARENG JAC BiCGStab double cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_0_q2("PARENG JAC BiCGStab double cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif

#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q1("PARENG JAC BiCGStab double sse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q1("PARENG JAC BiCGStab double sse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q1("PARENG JAC BiCGStab double sse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_1_q1("PARENG JAC BiCGStab double sse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_1_q2("PARENG JAC BiCGStab double sse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q2("PARENG JAC BiCGStab double sse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q2("PARENG JAC BiCGStab double sse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q2("PARENG JAC BiCGStab double sse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q1("PARENG JAC BiCGStab double mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q1("PARENG JAC BiCGStab double mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_1_q1("PARENG JAC BiCGStab double mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_1_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q2("PARENG JAC BiCGStab double mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q2("PARENG JAC BiCGStab double mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q2("PARENG JAC BiCGStab double mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q1("PARENG JAC BiCGStab double cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q1("PARENG JAC BiCGStab double cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_1_q1("PARENG JAC BiCGStab double cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_1_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_1_q2("PARENG JAC BiCGStab double cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_1_q2("PARENG JAC BiCGStab double cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_1_q2("PARENG JAC BiCGStab double cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif

#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q1("PARENG JAC BiCGStab double sse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q1("PARENG JAC BiCGStab double sse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q1("PARENG JAC BiCGStab double sse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_2_q1("PARENG JAC BiCGStab double sse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_2_q2("PARENG JAC BiCGStab double sse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q2("PARENG JAC BiCGStab double sse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q2("PARENG JAC BiCGStab double sse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q2("PARENG JAC BiCGStab double sse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q1("PARENG JAC BiCGStab double mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q1("PARENG JAC BiCGStab double mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_2_q1("PARENG JAC BiCGStab double mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_2_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q2("PARENG JAC BiCGStab double mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q2("PARENG JAC BiCGStab double mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q2("PARENG JAC BiCGStab double mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q1("PARENG JAC BiCGStab double cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q1("PARENG JAC BiCGStab double cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_2_q1("PARENG JAC BiCGStab double cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_2_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_2_q2("PARENG JAC BiCGStab double cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_2_q2("PARENG JAC BiCGStab double cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_2_q2("PARENG JAC BiCGStab double cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif


#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q1("PARENG JAC BiCGStab double sse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q1("PARENG JAC BiCGStab double sse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q1("PARENG JAC BiCGStab double sse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_3_q1("PARENG JAC BiCGStab double sse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_3_q2("PARENG JAC BiCGStab double sse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q2("PARENG JAC BiCGStab double sse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q2("PARENG JAC BiCGStab double sse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q2("PARENG JAC BiCGStab double sse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q1("PARENG JAC BiCGStab double mcsse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q1("PARENG JAC BiCGStab double mcsse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_3_q1("PARENG JAC BiCGStab double mcsse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_3_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q2("PARENG JAC BiCGStab double mcsse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q2("PARENG JAC BiCGStab double mcsse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q2("PARENG JAC BiCGStab double mcsse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q1("PARENG JAC BiCGStab double cuda L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q1("PARENG JAC BiCGStab double cuda L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_3_q1("PARENG JAC BiCGStab double cuda L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_3_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_3_q2("PARENG JAC BiCGStab double cuda L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_3_q2("PARENG JAC BiCGStab double cuda L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_3_q2("PARENG JAC BiCGStab double cuda L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif


#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q1("PARENG JAC BiCGStab double sse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q1("PARENG JAC BiCGStab double sse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q1("PARENG JAC BiCGStab double sse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_4_q1("PARENG JAC BiCGStab double sse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_4_q2("PARENG JAC BiCGStab double sse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q2("PARENG JAC BiCGStab double sse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q2("PARENG JAC BiCGStab double sse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::JAC> sse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q2("PARENG JAC BiCGStab double sse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q1("PARENG JAC BiCGStab double mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q1("PARENG JAC BiCGStab double mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_4_q1("PARENG JAC BiCGStab double mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_4_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q2("PARENG JAC BiCGStab double mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q2("PARENG JAC BiCGStab double mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q2("PARENG JAC BiCGStab double mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q1("PARENG JAC BiCGStab double cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q1("PARENG JAC BiCGStab double cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_10_4_q1("PARENG JAC BiCGStab double cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_6_4_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_7_4_q2("PARENG JAC BiCGStab double cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_8_4_q2("PARENG JAC BiCGStab double cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced_bench_bicgstabjac_sparse_prolmat_double_9_4_q2("PARENG JAC BiCGStab double cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif

template <typename Tag_, typename DT1_>
class BiCGStabPoissonAdvancedEllBENCH<Tag_, DT1_, methods::SPAI>:
    public Benchmark
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;

    public:
        BiCGStabPoissonAdvancedEllBENCH(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
        }

        virtual void run()
        {
            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";

            std::string A_file(file_base);
            A_file += "A_";
            A_file += stringify(_size);
            A_file += ".ell";
            SparseMatrixELL<DT1_> system(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

            std::string rhs_file(file_base);
            rhs_file += "rhs_";
            rhs_file += stringify(_size);
            DenseVector<DT1_> rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0));

            std::string precon_file(file_base);
            precon_file += "A_";
            precon_file += stringify(_size);
            precon_file += stringify("_spai");
            precon_file += ".ell";
            SparseMatrixELL<DT1_> precon = MatrixIO<io_formats::ELL>::read_matrix(precon_file, DT1_(0));

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            unsigned long used_iters(0);

            BENCHMARK(
                      (PBiCGStab<Tag_, methods::VAR>::value(system, rhs, result, precon, 5000ul, used_iters));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q1("PARENG SPAI BiCGStab double sse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q1("PARENG SPAI BiCGStab double sse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q1("PARENG SPAI BiCGStab double sse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_0_q1("PARENG SPAI BiCGStab double sse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_0_q2("PARENG SPAI BiCGStab double sse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q2("PARENG SPAI BiCGStab double sse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q2("PARENG SPAI BiCGStab double sse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q2("PARENG SPAI BiCGStab double sse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");


BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q1("PARENG SPAI BiCGStab double mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q1("PARENG SPAI BiCGStab double mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_0_q1("PARENG SPAI BiCGStab double mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_0_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q2("PARENG SPAI BiCGStab double mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q2("PARENG SPAI BiCGStab double mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q2("PARENG SPAI BiCGStab double mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q1("PARENG SPAI BiCGStab double cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q1("PARENG SPAI BiCGStab double cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_0_q1("PARENG SPAI BiCGStab double cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_0_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_0_q2("PARENG SPAI BiCGStab double cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_0_q2("PARENG SPAI BiCGStab double cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_0_q2("PARENG SPAI BiCGStab double cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif

#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q1("PARENG SPAI BiCGStab double sse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q1("PARENG SPAI BiCGStab double sse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q1("PARENG SPAI BiCGStab double sse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_1_q1("PARENG SPAI BiCGStab double sse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_1_q2("PARENG SPAI BiCGStab double sse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q2("PARENG SPAI BiCGStab double sse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q2("PARENG SPAI BiCGStab double sse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q2("PARENG SPAI BiCGStab double sse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q1("PARENG SPAI BiCGStab double mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q1("PARENG SPAI BiCGStab double mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_1_q1("PARENG SPAI BiCGStab double mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_1_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q2("PARENG SPAI BiCGStab double mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q2("PARENG SPAI BiCGStab double mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q2("PARENG SPAI BiCGStab double mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q1("PARENG SPAI BiCGStab double cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q1("PARENG SPAI BiCGStab double cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_1_q1("PARENG SPAI BiCGStab double cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_1_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_1_q2("PARENG SPAI BiCGStab double cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_1_q2("PARENG SPAI BiCGStab double cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_1_q2("PARENG SPAI BiCGStab double cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif

#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q1("PARENG SPAI BiCGStab double sse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q1("PARENG SPAI BiCGStab double sse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q1("PARENG SPAI BiCGStab double sse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_2_q1("PARENG SPAI BiCGStab double sse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_2_q2("PARENG SPAI BiCGStab double sse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q2("PARENG SPAI BiCGStab double sse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q2("PARENG SPAI BiCGStab double sse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q2("PARENG SPAI BiCGStab double sse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q1("PARENG SPAI BiCGStab double mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q1("PARENG SPAI BiCGStab double mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_2_q1("PARENG SPAI BiCGStab double mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_2_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q2("PARENG SPAI BiCGStab double mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q2("PARENG SPAI BiCGStab double mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q2("PARENG SPAI BiCGStab double mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q1("PARENG SPAI BiCGStab double cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q1("PARENG SPAI BiCGStab double cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_2_q1("PARENG SPAI BiCGStab double cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_2_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_2_q2("PARENG SPAI BiCGStab double cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_2_q2("PARENG SPAI BiCGStab double cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_2_q2("PARENG SPAI BiCGStab double cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif


#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q1("PARENG SPAI BiCGStab double sse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q1("PARENG SPAI BiCGStab double sse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q1("PARENG SPAI BiCGStab double sse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_3_q1("PARENG SPAI BiCGStab double sse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_3_q2("PARENG SPAI BiCGStab double sse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q2("PARENG SPAI BiCGStab double sse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q2("PARENG SPAI BiCGStab double sse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q2("PARENG SPAI BiCGStab double sse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q1("PARENG SPAI BiCGStab double mcsse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q1("PARENG SPAI BiCGStab double mcsse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_3_q1("PARENG SPAI BiCGStab double mcsse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_3_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q2("PARENG SPAI BiCGStab double mcsse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q2("PARENG SPAI BiCGStab double mcsse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q2("PARENG SPAI BiCGStab double mcsse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q1("PARENG SPAI BiCGStab double cuda L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q1("PARENG SPAI BiCGStab double cuda L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_3_q1("PARENG SPAI BiCGStab double cuda L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_3_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_3_q2("PARENG SPAI BiCGStab double cuda L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_3_q2("PARENG SPAI BiCGStab double cuda L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_3_q2("PARENG SPAI BiCGStab double cuda L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif


#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q1("PARENG SPAI BiCGStab double sse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q1("PARENG SPAI BiCGStab double sse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q1("PARENG SPAI BiCGStab double sse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_4_q1("PARENG SPAI BiCGStab double sse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_4_q2("PARENG SPAI BiCGStab double sse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q2("PARENG SPAI BiCGStab double sse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q2("PARENG SPAI BiCGStab double sse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::SSE, double, methods::SPAI> sse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q2("PARENG SPAI BiCGStab double sse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");

BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q1("PARENG SPAI BiCGStab double mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q1("PARENG SPAI BiCGStab double mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_4_q1("PARENG SPAI BiCGStab double mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_4_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q2("PARENG SPAI BiCGStab double mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q2("PARENG SPAI BiCGStab double mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q2("PARENG SPAI BiCGStab double mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q1("PARENG SPAI BiCGStab double cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q1("PARENG SPAI BiCGStab double cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_10_4_q1("PARENG SPAI BiCGStab double cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_6_4_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_7_4_q2("PARENG SPAI BiCGStab double cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_8_4_q2("PARENG SPAI BiCGStab double cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced_bench_bicgstab_sparse_prolmat_double_9_4_q2("PARENG SPAI BiCGStab double cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_");
#endif
#endif
