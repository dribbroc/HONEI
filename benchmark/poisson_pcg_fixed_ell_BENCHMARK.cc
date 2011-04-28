/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

//#define SOLVER_VERBOSE 1
#include <honei/math/conjugate_gradients.hh>
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

template <typename Tag_, typename DT1_>
class CGPoissonAdvancedEllBENCH:
    public Benchmark
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;
        unsigned _nc;

        static unsigned long _level_to_size(unsigned long level, unsigned elem_type) // 0 = q1 , 1 = q1t , 2 = q2
        {
            switch(level)
            {
                case 10:
                    {
                        if(elem_type == 0)
                            return 2101248;
                    }
                case 9:
                    {
                        if(elem_type == 0)
                            return 526336;
                        else if(elem_type == 1)
                            return 1050624;
                        else
                            return 2101248;
                    }
                case 8:
                    {
                        if(elem_type == 0)
                            return 132096;
                        else if(elem_type == 1)
                            return 263168;
                        else
                            return 526336;
                    }
                case 7:
                    {
                        if(elem_type == 0)
                            return 33280;
                        else if(elem_type == 1)
                            return 66048;
                        else
                            return 132096;
                    }
                case 6:
                    {
                        if(elem_type == 0)
                            return 8448;
                        else if(elem_type == 1)
                            return 16640;
                        else
                            return 33280;
                    }
                case 5:
                    {
                        if(elem_type == 0)
                            return 2176;
                        else if(elem_type == 1)
                            return 4224;
                        else
                            return 8448;
                    }
                case 4:
                    {
                        if(elem_type == 0)
                            return 576;
                        else if(elem_type == 1)
                            return 1088;
                        else
                            return 2176;
                    }
                case 3:
                    {
                        if(elem_type == 0)
                            return 160;
                        else if(elem_type == 1)
                            return 288;
                        else
                            return 576;
                    }
                case 2:
                    {
                        if(elem_type == 0)
                            return 48;
                        else if(elem_type == 1)
                            return 80;
                        else
                            return 160;
                    }
                case 1:
                    {
                        if(elem_type == 0)
                            return 16;
                        else if(elem_type == 1)
                            return 24;
                        else
                            return 48;
                    }
                default:
                    return 1;
            }
        }


    public:
        CGPoissonAdvancedEllBENCH(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
            _nc = nc;
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

            DenseVector<DT1_> diag(_level_to_size(_size, _nc));
            for(unsigned long j(0) ; j < _level_to_size(_size, _nc) ; ++j)
                diag[j] = system(j, j);

            ElementInverse<Tag_>::value(diag);
            Scale<Tag_>::value(diag, 0.7);

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            unsigned long used_iters(0);

            BENCHMARK(
                      (ConjugateGradients<Tag_, methods::JAC>::value(system, rhs, result, diag, 5000ul, used_iters));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q1("PARENG CG double sse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q1("PARENG CG double sse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q1("PARENG CG double sse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_0_q1("PARENG CG double sse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_6_0_q2("PARENG CG double sse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q2("PARENG CG double sse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q2("PARENG CG double sse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q2("PARENG CG double sse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q1t("PARENG CG double sse L7, q1t sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q1t("PARENG CG double sse L8, q1t sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q1t("PARENG CG double sse L9, q1t sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_0_q1t("PARENG CG double sse L10, q1t sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);



CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q1("PARENG CG double mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q1("PARENG CG double mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q1("PARENG CG double mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_10_0_q1("PARENG CG double mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_6_0_q2("PARENG CG double mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q2("PARENG CG double mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q2("PARENG CG double mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q2("PARENG CG double mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q1("PARENG CG double cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q1("PARENG CG double cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q1("PARENG CG double cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_10_0_q1("PARENG CG double cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_6_0_q2("PARENG CG double cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_0_q2("PARENG CG double cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_0_q2("PARENG CG double cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_0_q2("PARENG CG double cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif

#ifdef HONEI_SSE
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q1("PARENG CG double sse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q1("PARENG CG double sse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q1("PARENG CG double sse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_1_q1("PARENG CG double sse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_6_1_q2("PARENG CG double sse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q2("PARENG CG double sse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q2("PARENG CG double sse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q2("PARENG CG double sse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q1("PARENG CG double mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q1("PARENG CG double mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q1("PARENG CG double mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_10_1_q1("PARENG CG double mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_6_1_q2("PARENG CG double mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q2("PARENG CG double mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q2("PARENG CG double mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q2("PARENG CG double mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q1("PARENG CG double cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q1("PARENG CG double cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q1("PARENG CG double cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_10_1_q1("PARENG CG double cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_6_1_q2("PARENG CG double cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_1_q2("PARENG CG double cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_1_q2("PARENG CG double cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_1_q2("PARENG CG double cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif

#ifdef HONEI_SSE
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q1("PARENG CG double sse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q1("PARENG CG double sse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q1("PARENG CG double sse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_2_q1("PARENG CG double sse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_6_2_q2("PARENG CG double sse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q2("PARENG CG double sse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q2("PARENG CG double sse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q2("PARENG CG double sse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q1("PARENG CG double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q1("PARENG CG double mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q1("PARENG CG double mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_10_2_q1("PARENG CG double mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_6_2_q2("PARENG CG double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q2("PARENG CG double mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q2("PARENG CG double mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q2("PARENG CG double mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q1("PARENG CG double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q1("PARENG CG double cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q1("PARENG CG double cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_10_2_q1("PARENG CG double cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_6_2_q2("PARENG CG double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_2_q2("PARENG CG double cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_2_q2("PARENG CG double cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_2_q2("PARENG CG double cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif


#ifdef HONEI_SSE
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q1("PARENG CG double sse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q1("PARENG CG double sse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q1("PARENG CG double sse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_3_q1("PARENG CG double sse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_6_3_q2("PARENG CG double sse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q2("PARENG CG double sse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q2("PARENG CG double sse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q2("PARENG CG double sse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q1("PARENG CG double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q1("PARENG CG double mcsse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q1("PARENG CG double mcsse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_10_3_q1("PARENG CG double mcsse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_6_3_q2("PARENG CG double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q2("PARENG CG double mcsse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q2("PARENG CG double mcsse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q2("PARENG CG double mcsse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q1("PARENG CG double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q1("PARENG CG double cuda L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q1("PARENG CG double cuda L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_10_3_q1("PARENG CG double cuda L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_6_3_q2("PARENG CG double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_3_q2("PARENG CG double cuda L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_3_q2("PARENG CG double cuda L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_3_q2("PARENG CG double cuda L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif


#ifdef HONEI_SSE
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q1("PARENG CG double sse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q1("PARENG CG double sse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q1("PARENG CG double sse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_10_4_q1("PARENG CG double sse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_6_4_q2("PARENG CG double sse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q2("PARENG CG double sse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q2("PARENG CG double sse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q2("PARENG CG double sse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q1("PARENG CG double mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q1("PARENG CG double mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q1("PARENG CG double mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_10_4_q1("PARENG CG double mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_6_4_q2("PARENG CG double mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q2("PARENG CG double mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q2("PARENG CG double mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q2("PARENG CG double mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q1("PARENG CG double cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q1("PARENG CG double cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q1("PARENG CG double cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_10_4_q1("PARENG CG double cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_6_4_q2("PARENG CG double cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_7_4_q2("PARENG CG double cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_8_4_q2("PARENG CG double cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedEllBENCH<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_sparse_prolmat_double_9_4_q2("PARENG CG double cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif

template <typename Tag_, typename DT1_>
class CGPoissonAdvancedELLBENCHSPAI:
    public Benchmark
{
    private:
        unsigned long _size;
        std::string _res_f, _file_base;
        unsigned long _sorting;
        unsigned _nc;

        static unsigned long _level_to_size(unsigned long level, unsigned elem_type) // 0 = q1 , 1 = q1t , 2 = q2
        {
            switch(level)
            {
                case 10:
                    {
                        if(elem_type == 0)
                            return 2101248;
                    }
                case 9:
                    {
                        if(elem_type == 0)
                            return 526336;
                        else if(elem_type == 1)
                            return 1050624;
                        else
                            return 2101248;
                    }
                case 8:
                    {
                        if(elem_type == 0)
                            return 132096;
                        else if(elem_type == 1)
                            return 263168;
                        else
                            return 526336;
                    }
                case 7:
                    {
                        if(elem_type == 0)
                            return 33280;
                        else if(elem_type == 1)
                            return 66048;
                        else
                            return 132096;
                    }
                case 6:
                    {
                        if(elem_type == 0)
                            return 8448;
                        else if(elem_type == 1)
                            return 16640;
                        else
                            return 33280;
                    }
                case 5:
                    {
                        if(elem_type == 0)
                            return 2176;
                        else if(elem_type == 1)
                            return 4224;
                        else
                            return 8448;
                    }
                case 4:
                    {
                        if(elem_type == 0)
                            return 576;
                        else if(elem_type == 1)
                            return 1088;
                        else
                            return 2176;
                    }
                case 3:
                    {
                        if(elem_type == 0)
                            return 160;
                        else if(elem_type == 1)
                            return 288;
                        else
                            return 576;
                    }
                case 2:
                    {
                        if(elem_type == 0)
                            return 48;
                        else if(elem_type == 1)
                            return 80;
                        else
                            return 160;
                    }
                case 1:
                    {
                        if(elem_type == 0)
                            return 16;
                        else if(elem_type == 1)
                            return 24;
                        else
                            return 48;
                    }
                default:
                    return 1;
            }
        }


    public:
        CGPoissonAdvancedELLBENCHSPAI(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
            _nc = nc;
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

            std::string spai_file(file_base);
            spai_file += "A_";
            spai_file += stringify(_size);
            spai_file += "_spai.ell";
            SparseMatrixELL<DT1_> precon(MatrixIO<io_formats::ELL>::read_matrix(spai_file, DT1_(0)));

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            unsigned long used_iters(0);

            BENCHMARK(
                      (ConjugateGradients<Tag_, methods::SPAI>::value(system, rhs, result, precon, 5000ul, used_iters));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q1("PARENG CG SPAI double sse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q1("PARENG CG SPAI double sse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q1("PARENG CG SPAI double sse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_0_q1("PARENG CG SPAI double sse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_0_q2("PARENG CG SPAI double sse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q2("PARENG CG SPAI double sse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q2("PARENG CG SPAI double sse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q2("PARENG CG SPAI double sse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q1t("PARENG CG SPAI double sse L7, q1t sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q1t("PARENG CG SPAI double sse L8, q1t sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q1t("PARENG CG SPAI double sse L9, q1t sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_0_q1t("PARENG CG SPAI double sse L10, q1t sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/q1t_sort_", 1);



CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q1("PARENG CG SPAI double mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q1("PARENG CG SPAI double mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q1("PARENG CG SPAI double mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_0_q1("PARENG CG SPAI double mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_0_q2("PARENG CG SPAI double mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q2("PARENG CG SPAI double mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q2("PARENG CG SPAI double mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q2("PARENG CG SPAI double mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q1("PARENG CG SPAI double cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q1("PARENG CG SPAI double cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q1("PARENG CG SPAI double cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_0_q1("PARENG CG SPAI double cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_0_q2("PARENG CG SPAI double cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_0_q2("PARENG CG SPAI double cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_0_q2("PARENG CG SPAI double cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_0_q2("PARENG CG SPAI double cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif

#ifdef HONEI_SSE
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q1("PARENG CG SPAI double sse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q1("PARENG CG SPAI double sse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q1("PARENG CG SPAI double sse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_1_q1("PARENG CG SPAI double sse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_1_q2("PARENG CG SPAI double sse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q2("PARENG CG SPAI double sse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q2("PARENG CG SPAI double sse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q2("PARENG CG SPAI double sse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q1("PARENG CG SPAI double mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q1("PARENG CG SPAI double mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q1("PARENG CG SPAI double mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_1_q1("PARENG CG SPAI double mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_1_q2("PARENG CG SPAI double mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q2("PARENG CG SPAI double mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q2("PARENG CG SPAI double mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q2("PARENG CG SPAI double mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q1("PARENG CG SPAI double cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q1("PARENG CG SPAI double cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q1("PARENG CG SPAI double cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_1_q1("PARENG CG SPAI double cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_1_q2("PARENG CG SPAI double cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_1_q2("PARENG CG SPAI double cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_1_q2("PARENG CG SPAI double cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_1_q2("PARENG CG SPAI double cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif

#ifdef HONEI_SSE
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q1("PARENG CG SPAI double sse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q1("PARENG CG SPAI double sse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q1("PARENG CG SPAI double sse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_2_q1("PARENG CG SPAI double sse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_2_q2("PARENG CG SPAI double sse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q2("PARENG CG SPAI double sse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q2("PARENG CG SPAI double sse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q2("PARENG CG SPAI double sse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q1("PARENG CG SPAI double mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q1("PARENG CG SPAI double mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q1("PARENG CG SPAI double mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_2_q1("PARENG CG SPAI double mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_2_q2("PARENG CG SPAI double mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q2("PARENG CG SPAI double mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q2("PARENG CG SPAI double mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q2("PARENG CG SPAI double mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q1("PARENG CG SPAI double cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q1("PARENG CG SPAI double cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q1("PARENG CG SPAI double cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_2_q1("PARENG CG SPAI double cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_2_q2("PARENG CG SPAI double cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_2_q2("PARENG CG SPAI double cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_2_q2("PARENG CG SPAI double cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_2_q2("PARENG CG SPAI double cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif


#ifdef HONEI_SSE
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q1("PARENG CG SPAI double sse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q1("PARENG CG SPAI double sse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q1("PARENG CG SPAI double sse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_3_q1("PARENG CG SPAI double sse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_3_q2("PARENG CG SPAI double sse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q2("PARENG CG SPAI double sse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q2("PARENG CG SPAI double sse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q2("PARENG CG SPAI double sse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q1("PARENG CG SPAI double mcsse L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q1("PARENG CG SPAI double mcsse L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q1("PARENG CG SPAI double mcsse L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_3_q1("PARENG CG SPAI double mcsse L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_3_q2("PARENG CG SPAI double mcsse L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q2("PARENG CG SPAI double mcsse L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q2("PARENG CG SPAI double mcsse L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q2("PARENG CG SPAI double mcsse L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q1("PARENG CG SPAI double cuda L7, q1 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q1("PARENG CG SPAI double cuda L8, q1 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q1("PARENG CG SPAI double cuda L9, q1 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_3_q1("PARENG CG SPAI double cuda L10, q1 sort 3", 10ul, 3ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_3_q2("PARENG CG SPAI double cuda L6, q2 sort 3", 6ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_3_q2("PARENG CG SPAI double cuda L7, q2 sort 3", 7ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_3_q2("PARENG CG SPAI double cuda L8, q2 sort 3", 8ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_3_q2("PARENG CG SPAI double cuda L9, q2 sort 3", 9ul, 3ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif


#ifdef HONEI_SSE
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q1("PARENG CG SPAI double sse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q1("PARENG CG SPAI double sse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q1("PARENG CG SPAI double sse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_4_q1("PARENG CG SPAI double sse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_4_q2("PARENG CG SPAI double sse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q2("PARENG CG SPAI double sse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q2("PARENG CG SPAI double sse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::SSE, double> sse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q2("PARENG CG SPAI double sse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q1("PARENG CG SPAI double mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q1("PARENG CG SPAI double mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q1("PARENG CG SPAI double mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_4_q1("PARENG CG SPAI double mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_4_q2("PARENG CG SPAI double mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q2("PARENG CG SPAI double mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q2("PARENG CG SPAI double mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::CPU::MultiCore::SSE, double> mcsse_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q2("PARENG CG SPAI double mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q1("PARENG CG SPAI double cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q1("PARENG CG SPAI double cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q1("PARENG CG SPAI double cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_10_4_q1("PARENG CG SPAI double cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_6_4_q2("PARENG CG SPAI double cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_7_4_q2("PARENG CG SPAI double cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_8_4_q2("PARENG CG SPAI double cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
CGPoissonAdvancedELLBENCHSPAI<tags::GPU::CUDA, double> cuda_poisson_advanced_bench_cg_spai_sparse_prolmat_double_9_4_q2("PARENG CG SPAI double cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif
