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
#include <honei/math/multigrid.hh>
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

template <typename ITag_, typename OTag_, typename DT1_>
class PoissonAdvancedBENCHMGSparseELLProlMat:
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
        PoissonAdvancedBENCHMGSparseELLProlMat(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            Benchmark(tag)
        {
            register_tag(OTag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
            _nc = nc;
        }

        virtual void run()
        {
            //unsigned long n(_level_to_size(_size, _nc));
            MGInfo<float, SparseMatrixELL<float> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;

            for(int i(0); i < 8; ++i)
            {
                (info.macro_border_mask)[i] = 2;
            }
            //set Neumann boundaries:
            //(info.macro_border_mask)[5] =1;


            info.min_level = 1;
            info.max_level = _size;
            info.n_max_iter = 3;
            info.initial_zero = true;
            info.tolerance = 1e-7;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            //info.n_max_iter_coarse = ((unsigned long)sqrt((float)(pow((float)2 , (float)info.max_level) + 1)*(pow((float)2 , (float)info.max_level) + 1)));
            info.n_max_iter_coarse = 1000;
            info.tolerance_coarse = 1e-2;
            info.adapt_correction_factor = 1.;

            std::cout << "Info is set up." << std::endl;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(i == 0)
                    size = 9;

                DenseVector<float> dummy_band(size, float(0));
                BandedMatrixQx<Q1Type, float> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<float> sm(ac_a);
                SparseMatrixELL<float> ac_s(sm);
                info.a.push_back(ac_s);
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_rhs(size, float(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                info.diags_inverted.push_back(dummy_band.copy());
            }

            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                std::cout << size << std::endl;
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                DenseVector<float> dummy_band(size, float(0));
                //info.diags_inverted.push_back(dummy_band.copy());
            }

            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";
            std::cout << "File:" << file_base << std::endl;
            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N(_level_to_size(i, _nc));
                DenseVector<float> current_rhs(N);
                std::string A_file(file_base);
                A_file += "A_";
                A_file += stringify(i);
                A_file += ".ell";
                SparseMatrixELL<float> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, float(0)));

                std::string rhs_file(file_base);
                rhs_file += "rhs_" + stringify(_size);
                if(i == info.max_level)
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, float(0));

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<float> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<OTag_>::value(scaled_diag_inverted);
                Scale<OTag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<float> prol(1, 1);
                        SparseMatrixELL<float> prolmat(prol);
                        info.prolmats.push_back(prolmat);
                        info.resmats.push_back(prolmat);
                    }
                    else
                    {
                        std::string prol_file(file_base);
                        prol_file += "prol_";
                        prol_file += stringify(i);
                        prol_file += ".ell";
                        SparseMatrixELL<float> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, float(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<float> prol(prolmat);
                        SparseMatrix<float> res(prol.columns(), prol.rows());
                        Transposition<OTag_>::value(prol, res);
                        SparseMatrixELL<float> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
                info.x[i] = null.copy();
            }

            /*for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                std::cout << "VECSIZE d " << info.d.at(i).size() << std::endl;
                std::cout << "VECSIZE rhs " << info.rhs.at(i).size() << std::endl;
                std::cout <<"SYSTEM"<<std::endl;
                std::cout <<"ROW: " << info.a.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.a.at(i).columns() << std::endl;
                std::cout <<"RESTRICTION"<<std::endl;
                std::cout <<"ROW: " << info.resmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.resmats.at(i).columns() << std::endl;
                std::cout <<"PROLONGATION"<<std::endl;
                std::cout <<"ROW: " << info.prolmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.prolmats.at(i).columns() << std::endl;

            }*/

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            std::string rhs_file(file_base);
            rhs_file += "rhs_" + stringify(_size);
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0)));

            std::string A_file(file_base);
            A_file += "A_";
            A_file += stringify(_size);
            A_file += ".ell";
            SparseMatrixELL<DT1_> system(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

            std::cout << "Outer data is set up." << std::endl;
            std::cout.flush();

            BENCHMARK(
                      (Multigrid<ITag_, OTag_, methods::PROLMAT, methods::JAC, methods::CYCLE::V, methods::MIXED >::value(system, rhs, result, (unsigned long)11, 1e-8, info));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE

//SORT 0
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_0_q1("PARENG MG double cuda-mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_0_q1("PARENG MG double cuda-mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_0_q1("PARENG MG double cuda-mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_10_0_q1("PARENG MG double cuda-mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_6_0_q2("PARENG MG double cuda-mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_0_q2("PARENG MG double cuda-mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_0_q2("PARENG MG double cuda-mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_0_q2("PARENG MG double cuda-mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_0_q1("PARENG MG double cuda-cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_0_q1("PARENG MG double cuda-cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_0_q1("PARENG MG double cuda-cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_10_0_q1("PARENG MG double cuda-cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_6_0_q2("PARENG MG double cuda-cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_0_q2("PARENG MG double cuda-cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_0_q2("PARENG MG double cuda-cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_0_q2("PARENG MG double cuda-cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 1
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_1_q1("PARENG MG double cuda-mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_1_q1("PARENG MG double cuda-mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_1_q1("PARENG MG double cuda-mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_10_1_q1("PARENG MG double cuda-mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_6_1_q2("PARENG MG double cuda-mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_1_q2("PARENG MG double cuda-mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_1_q2("PARENG MG double cuda-mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_1_q2("PARENG MG double cuda-mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_1_q1("PARENG MG double cuda-cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_1_q1("PARENG MG double cuda-cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_1_q1("PARENG MG double cuda-cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_10_1_q1("PARENG MG double cuda-cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_6_1_q2("PARENG MG double cuda-cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_1_q2("PARENG MG double cuda-cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_1_q2("PARENG MG double cuda-cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_1_q2("PARENG MG double cuda-cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 2
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_2_q1("PARENG MG double cuda-mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_2_q1("PARENG MG double cuda-mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_2_q1("PARENG MG double cuda-mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_10_2_q1("PARENG MG double cuda-mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_6_2_q2("PARENG MG double cuda-mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_2_q2("PARENG MG double cuda-mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_2_q2("PARENG MG double cuda-mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_2_q2("PARENG MG double cuda-mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_2_q1("PARENG MG double cuda-cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_2_q1("PARENG MG double cuda-cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_2_q1("PARENG MG double cuda-cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_10_2_q1("PARENG MG double cuda-cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_6_2_q2("PARENG MG double cuda-cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_2_q2("PARENG MG double cuda-cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_2_q2("PARENG MG double cuda-cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_2_q2("PARENG MG double cuda-cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 4
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_4_q1("PARENG MG double cuda-mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_4_q1("PARENG MG double cuda-mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_4_q1("PARENG MG double cuda-mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_10_4_q1("PARENG MG double cuda-mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_6_4_q2("PARENG MG double cuda-mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_7_4_q2("PARENG MG double cuda-mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_8_4_q2("PARENG MG double cuda-mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_double_9_4_q2("PARENG MG double cuda-mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_4_q1("PARENG MG double cuda-cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_4_q1("PARENG MG double cuda-cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_4_q1("PARENG MG double cuda-cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_10_4_q1("PARENG MG double cuda-cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_6_4_q2("PARENG MG double cuda-cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_7_4_q2("PARENG MG double cuda-cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_8_4_q2("PARENG MG double cuda-cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMat<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_double_9_4_q2("PARENG MG double cuda-cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif
#endif

//-----------------SPAI-----------------------------------

template <typename ITag_, typename OTag_, typename DT1_>
class PoissonAdvancedBENCHMGSparseELLProlMatSPAI:
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
        PoissonAdvancedBENCHMGSparseELLProlMatSPAI(const std::string & tag,
                unsigned long level, unsigned long sorting, std::string file_base, unsigned nc) :
            Benchmark(tag)
        {
            register_tag(OTag_::name);
            _size = level;
            _sorting = sorting;
            _file_base = file_base;
            _nc = nc;
        }

        virtual void run()
        {
            //unsigned long n(_level_to_size(_size, _nc));
            MGInfo<float, SparseMatrixELL<float> > info;
            //configuration constants: /TODO: set/allocate!!!
            info.is_smoother = false;

            for(int i(0); i < 8; ++i)
            {
                (info.macro_border_mask)[i] = 2;
            }
            //set Neumann boundaries:
            //(info.macro_border_mask)[5] =1;


            info.min_level = 1;
            info.max_level = _size;
            info.n_max_iter = 2;
            info.initial_zero = true;
            info.tolerance = 1e-8;
            info.convergence_check = true;

            info.n_pre_smooth = 4;
            info.n_post_smooth = 4;
            //info.n_max_iter_coarse = ((unsigned long)sqrt((float)(pow((float)2 , (float)info.max_level) + 1)*(pow((float)2 , (float)info.max_level) + 1)));
            info.n_max_iter_coarse = 10000;
            info.tolerance_coarse = 1e-8;
            info.adapt_correction_factor = 1.;

            std::cout << "Info is set up." << std::endl;

            for (unsigned long i(0) ; i < info.min_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(i == 0)
                    size = 9;

                DenseVector<float> dummy_band(size, float(0));
                BandedMatrixQx<Q1Type, float> ac_a(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                SparseMatrix<float> sm(ac_a);
                SparseMatrixELL<float> ac_s(sm);
                info.a.push_back(ac_s);
                info.prolmats.push_back(ac_s.copy());
                info.resmats.push_back(ac_s.copy());
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_rhs(size, float(0));
                info.rhs.push_back(ac_rhs);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);
                info.spais.push_back(ac_s.copy());
                info.temp_0.push_back(dummy_band.copy());
                info.temp_1.push_back(dummy_band.copy());

                info.diags_inverted.push_back(dummy_band.copy());
            }

            for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                std::cout << size << std::endl;
                // iteration vectors
                DenseVector<float> ac_c(size, float(0));
                info.c.push_back(ac_c);
                DenseVector<float> ac_d(size, float(0));
                info.d.push_back(ac_d);
                DenseVector<float> ac_x(size, float(0));
                info.x.push_back(ac_x);

                DenseVector<float> dummy_band(size, float(0));
                //info.diags_inverted.push_back(dummy_band.copy());
            }

            std::string file_base(HONEI_SOURCEDIR);
            file_base += _file_base + stringify(_sorting) + "/";
            std::cout << "File:" << file_base << std::endl;
            //assemble all needed levels' matrices:
            for(unsigned long i(info.min_level); i <= info.max_level; ++i)
            {
                unsigned long N(_level_to_size(i, _nc));
                DenseVector<float> current_rhs(N);
                std::string A_file(file_base);
                A_file += "A_";
                A_file += stringify(i);
                A_file += ".ell";
                SparseMatrixELL<float> smell(MatrixIO<io_formats::ELL>::read_matrix(A_file, float(0)));

                std::string rhs_file(file_base);
                rhs_file += "rhs_" + stringify(_size);
                if(i == info.max_level)
                    current_rhs = VectorIO<io_formats::EXP>::read_vector(rhs_file, float(0));

                info.rhs.push_back(current_rhs);
                info.a.push_back(smell);

                DenseVector<float> scaled_diag_inverted(N);
                for(unsigned long j(0) ; j < N ; ++ j)
                    scaled_diag_inverted[j] = smell(j, j);

                ElementInverse<OTag_>::value(scaled_diag_inverted);
                Scale<OTag_>::value(scaled_diag_inverted, 0.7);

                info.diags_inverted.push_back(scaled_diag_inverted.copy());

                info.temp_0.push_back(scaled_diag_inverted.copy());
                info.temp_1.push_back(scaled_diag_inverted.copy());
                SparseMatrix<float> spai(smell);
                std::string SPAI_file(file_base);
                SPAI_file += "A_";
                SPAI_file += stringify(i);
                SPAI_file += "_spai.ell";
                SparseMatrixELL<float> spai_m(MatrixIO<io_formats::ELL>::read_matrix(SPAI_file, float(0)));
                info.spais.push_back(spai_m);

                if(i >= info.min_level)
                {
                    if(i == 1)
                    {
                        SparseMatrix<float> prol(1, 1);
                        SparseMatrixELL<float> prolmat(prol);
                        info.prolmats.push_back(prolmat);
                        info.resmats.push_back(prolmat);
                    }
                    else
                    {
                        std::string prol_file(file_base);
                        prol_file += "prol_";
                        prol_file += stringify(i);
                        prol_file += ".ell";
                        SparseMatrixELL<float> prolmat(MatrixIO<io_formats::ELL>::read_matrix(prol_file, float(0)));
                        info.prolmats.push_back(prolmat);

                        SparseMatrix<float> prol(prolmat);
                        SparseMatrix<float> res(prol.columns(), prol.rows());
                        Transposition<OTag_>::value(prol, res);
                        SparseMatrixELL<float> resmat(res);
                        info.resmats.push_back(resmat);
                    }
                }
            }
            //clear x data
            for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                unsigned long size(_level_to_size(i, _nc));
                if(size==0)
                    size = 9;

                DenseVector<float> null(size , float(0));
                info.x[i] = null.copy();
            }

            /*for(unsigned long i(0) ; i < info.max_level ; ++i)
            {
                std::cout << "VECSIZE d " << info.d.at(i).size() << std::endl;
                std::cout << "VECSIZE rhs " << info.rhs.at(i).size() << std::endl;
                std::cout <<"SYSTEM"<<std::endl;
                std::cout <<"ROW: " << info.a.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.a.at(i).columns() << std::endl;
                std::cout <<"RESTRICTION"<<std::endl;
                std::cout <<"ROW: " << info.resmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.resmats.at(i).columns() << std::endl;
                std::cout <<"PROLONGATION"<<std::endl;
                std::cout <<"ROW: " << info.prolmats.at(i).rows() << std::endl;
                std::cout <<"COLS: " << info.prolmats.at(i).columns() << std::endl;

            }*/

            std::string init_file(file_base);
            init_file += "init_" + stringify(_size);
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(init_file, DT1_(0)));

            std::string rhs_file(file_base);
            rhs_file += "rhs_" + stringify(_size);
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(rhs_file, DT1_(0)));

            std::string A_file(file_base);
            A_file += "A_";
            A_file += stringify(_size);
            A_file += ".ell";
            SparseMatrixELL<DT1_> system(MatrixIO<io_formats::ELL>::read_matrix(A_file, DT1_(0)));

            std::cout << "Outer data is set up." << std::endl;
            std::cout.flush();

            BENCHMARK(
                      (Multigrid<ITag_, OTag_, methods::PROLMAT, methods::SPAI, methods::CYCLE::V, methods::MIXED >::value(system, rhs, result, (unsigned long)11, 1e-8, info));
                     );

            evaluate();
        }
};
#ifdef HONEI_SSE
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
//SORT 0
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_0_q1("PARENG MG SPAI double cuda-mcsse L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_0_q1("PARENG MG SPAI double cuda-mcsse L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_0_q1("PARENG MG SPAI double cuda-mcsse L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_0_q1("PARENG MG SPAI double cuda-mcsse L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_0_q2("PARENG MG SPAI double cuda-mcsse L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_0_q2("PARENG MG SPAI double cuda-mcsse L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_0_q2("PARENG MG SPAI double cuda-mcsse L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_0_q2("PARENG MG SPAI double cuda-mcsse L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_0_q1("PARENG MG SPAI double cuda-cuda L7, q1 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_0_q1("PARENG MG SPAI double cuda-cuda L8, q1 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_0_q1("PARENG MG SPAI double cuda-cuda L9, q1 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_0_q1("PARENG MG SPAI double cuda-cuda L10, q1 sort 0", 10ul, 0ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_0_q2("PARENG MG SPAI double cuda-cuda L6, q2 sort 0", 6ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_0_q2("PARENG MG SPAI double cuda-cuda L7, q2 sort 0", 7ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_0_q2("PARENG MG SPAI double cuda-cuda L8, q2 sort 0", 8ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_0_q2("PARENG MG SPAI double cuda-cuda L9, q2 sort 0", 9ul, 0ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 1
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_1_q1("PARENG MG SPAI double cuda-mcsse L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_1_q1("PARENG MG SPAI double cuda-mcsse L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_1_q1("PARENG MG SPAI double cuda-mcsse L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_1_q1("PARENG MG SPAI double cuda-mcsse L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_1_q2("PARENG MG SPAI double cuda-mcsse L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_1_q2("PARENG MG SPAI double cuda-mcsse L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_1_q2("PARENG MG SPAI double cuda-mcsse L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_1_q2("PARENG MG SPAI double cuda-mcsse L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_1_q1("PARENG MG SPAI double cuda-cuda L7, q1 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_1_q1("PARENG MG SPAI double cuda-cuda L8, q1 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_1_q1("PARENG MG SPAI double cuda-cuda L9, q1 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_1_q1("PARENG MG SPAI double cuda-cuda L10, q1 sort 1", 10ul, 1ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_1_q2("PARENG MG SPAI double cuda-cuda L6, q2 sort 1", 6ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_1_q2("PARENG MG SPAI double cuda-cuda L7, q2 sort 1", 7ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_1_q2("PARENG MG SPAI double cuda-cuda L8, q2 sort 1", 8ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_1_q2("PARENG MG SPAI double cuda-cuda L9, q2 sort 1", 9ul, 1ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 2
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_2_q1("PARENG MG SPAI double cuda-mcsse L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_2_q1("PARENG MG SPAI double cuda-mcsse L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_2_q1("PARENG MG SPAI double cuda-mcsse L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_2_q1("PARENG MG SPAI double cuda-mcsse L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_2_q2("PARENG MG SPAI double cuda-mcsse L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_2_q2("PARENG MG SPAI double cuda-mcsse L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_2_q2("PARENG MG SPAI double cuda-mcsse L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_2_q2("PARENG MG SPAI double cuda-mcsse L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_2_q1("PARENG MG SPAI double cuda-cuda L7, q1 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_2_q1("PARENG MG SPAI double cuda-cuda L8, q1 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_2_q1("PARENG MG SPAI double cuda-cuda L9, q1 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_2_q1("PARENG MG SPAI double cuda-cuda L10, q1 sort 2", 10ul, 2ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_2_q2("PARENG MG SPAI double cuda-cuda L6, q2 sort 2", 6ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_2_q2("PARENG MG SPAI double cuda-cuda L7, q2 sort 2", 7ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_2_q2("PARENG MG SPAI double cuda-cuda L8, q2 sort 2", 8ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_2_q2("PARENG MG SPAI double cuda-cuda L9, q2 sort 2", 9ul, 2ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);

//SORT 4
//CUDA-MCSSE
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_4_q1("PARENG MG SPAI double cuda-mcsse L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_4_q1("PARENG MG SPAI double cuda-mcsse L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_4_q1("PARENG MG SPAI double cuda-mcsse L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_4_q1("PARENG MG SPAI double cuda-mcsse L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_4_q2("PARENG MG SPAI double cuda-mcsse L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_4_q2("PARENG MG SPAI double cuda-mcsse L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_4_q2("PARENG MG SPAI double cuda-mcsse L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::CPU::MultiCore::SSE, double> cudamc_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_4_q2("PARENG MG SPAI double cuda-mcsse L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
//CUDA-CUDA
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_4_q1("PARENG MG SPAI double cuda-cuda L7, q1 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_4_q1("PARENG MG SPAI double cuda-cuda L8, q1 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_4_q1("PARENG MG SPAI double cuda-cuda L9, q1 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_10_4_q1("PARENG MG SPAI double cuda-cuda L10, q1 sort 4", 10ul, 4ul, "/honei/math/testdata/poisson_advanced/sort_", 0);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_6_4_q2("PARENG MG SPAI double cuda-cuda L6, q2 sort 4", 6ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_7_4_q2("PARENG MG SPAI double cuda-cuda L7, q2 sort 4", 7ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_8_4_q2("PARENG MG SPAI double cuda-cuda L8, q2 sort 4", 8ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
PoissonAdvancedBENCHMGSparseELLProlMatSPAI<tags::GPU::CUDA, tags::GPU::CUDA, double> cuda_poisson_advanced_bench_mg_sparse_prolmat_spai_double_9_4_q2("PARENG MG SPAI double cuda-cuda L9, q2 sort 4", 9ul, 4ul, "/honei/math/testdata/poisson_advanced/q2_sort_", 2);
#endif
#endif
#endif
