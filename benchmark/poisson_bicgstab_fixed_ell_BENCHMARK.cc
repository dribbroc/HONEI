/* vim: set sw=4 sts=4 et foldmethod=syntax : */

//#define SOLVER_VERBOSE 1
#include <honei/math/bicgstab.hh>
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
                tsmatrix3(i , i) = 1/system(i, i);
            }
            SparseMatrixELL<DT1_> smatrix3(tsmatrix3);

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            unsigned long used_iters(0);

            BENCHMARK(
                      (BiCGStabSolver<Tag_, methods::VAR>::value(system, smatrix3, rhs, result, 10000ul, used_iters, DT1_(1e-8)));
                     );

            std::cout << "Iters: " << used_iters << std::endl;

            evaluate();
        }
};
#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_2_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_3_q1("PARENG JAC BiCGStab double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_2_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::JAC> mcsse_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_3_q2("PARENG JAC BiCGStab double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_2_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_3_q1("PARENG JAC BiCGStab double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_2_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::JAC> cuda_poisson_advanced4_bench_bicgstabjac_sparse_prolmat_double_7_3_q2("PARENG JAC BiCGStab double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
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
                      (BiCGStabSolver<Tag_, methods::VAR>::value(system, precon, rhs, result, 10000ul, used_iters, DT1_(1e-8)));
                     );

            std::cout << "Iters: " << used_iters << std::endl;

            evaluate();
        }
};
#ifdef HONEI_SSE
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_7_2_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_7_3_q1("PARENG SPAI BiCGStab double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_6_2_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double, methods::SPAI> mcsse_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_6_3_q2("PARENG SPAI BiCGStab double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_7_2_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_7_3_q1("PARENG SPAI BiCGStab double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced4/sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_6_2_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
BiCGStabPoissonAdvancedEllBENCH<tags::GPU::CUDA, double, methods::SPAI> cuda_poisson_advanced4_bench_bicgstabspai_sparse_prolmat_double_6_3_q2("PARENG SPAI BiCGStab double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced4/q2_sort_");
#endif
#endif

