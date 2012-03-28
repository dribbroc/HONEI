/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmmund.de>
 * Copyright (c) 2012 Markus Geveler <admin@honei.org>
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
#include <honei/la/product.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/transposition.hh>
#include <benchmark/benchmark.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/backends/cuda/gpu_pool.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace std;


template <typename Tag_, typename DT_>
class Pareng:
    public Benchmark
{
    private:
        std::string _file_name;
        unsigned long _count;

    public:
        Pareng(const std::string & tag, std::string filename, unsigned long count) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _file_name = filename;
            _count = count;
        }

        virtual void run()
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/";
            _file_name = filebase + _file_name;
            SparseMatrixELL<DT_> smatrix = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            DenseVector<DT_> x(smatrix.columns());
            DenseVector<DT_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 100 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            {
            BenchmarkInfo info;
            info.flops = smatrix.used_elements() * 2;
            info.load = smatrix.used_elements() * 3 * sizeof(DT_);
            info.store = smatrix.used_elements() * 1 * sizeof(DT_);
            evaluate(info * 100);
            }
            std::cout<<"Non Zero Elements: "<<smatrix.used_elements()<<std::endl;
        }
};

#ifdef HONEI_SSE
#ifdef HONEI_CUDA_DOUBLE
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_0_7_double("mcsse q1 sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_0_8_double("mcsse q1 sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_0_9_double("mcsse q1 sort_0 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_1_7_double("mcsse q1 sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_1_8_double("mcsse q1 sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_1_9_double("mcsse q1 sort_1 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_2_7_double("mcsse q1 sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_2_8_double("mcsse q1 sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_2_9_double("mcsse q1 sort_2 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_3_7_double("mcsse q1 sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_3_8_double("mcsse q1 sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_3_9_double("mcsse q1 sort_3 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_4_7_double("mcsse q1 sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_4_8_double("mcsse q1 sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1_sort_4_9_double("mcsse q1 sort_4 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_0_6_double("mcsse q2 sort_0 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_0_7_double("mcsse q2 sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_0_8_double("mcsse q2 sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_1_6_double("mcsse q2 sort_1 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_1_7_double("mcsse q2 sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_1_8_double("mcsse q2 sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_2_6_double("mcsse q2 sort_2 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_2_7_double("mcsse q2 sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_2_8_double("mcsse q2 sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_3_6_double("mcsse q2 sort_3 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_3_7_double("mcsse q2 sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_3_8_double("mcsse q2 sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_4_6_double("mcsse q2 sort_4 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_4_7_double("mcsse q2 sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q2_sort_4_8_double("mcsse q2 sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_0_6_double("mcsse q1t sort_0 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_0_7_double("mcsse q1t sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_0_8_double("mcsse q1t sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_1_6_double("mcsse q1t sort_1 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_1_7_double("mcsse q1t sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_1_8_double("mcsse q1t sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_2_6_double("mcsse q1t sort_2 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_2_7_double("mcsse q1t sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_2_8_double("mcsse q1t sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_3_6_double("mcsse q1t sort_3 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_3_7_double("mcsse q1t sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_3_8_double("mcsse q1t sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_4_6_double("mcsse q1t sort_4 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_4_7_double("mcsse q1t sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_q1t_sort_4_8_double("mcsse q1t sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_0_7_double("cuda q1 sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_0_8_double("cuda q1 sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_0_9_double("cuda q1 sort_0 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_0/A_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_1_7_double("cuda q1 sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_1_8_double("cuda q1 sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_1_9_double("cuda q1 sort_1 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_1/A_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_2_7_double("cuda q1 sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_2_8_double("cuda q1 sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_2_9_double("cuda q1 sort_2 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_2/A_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_3_7_double("cuda q1 sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_3_8_double("cuda q1 sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_3_9_double("cuda q1 sort_3 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_3/A_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_4_7_double("cuda q1 sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_4_8_double("cuda q1 sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1_sort_4_9_double("cuda q1 sort_4 double L9 ELLPACK stiffness", "testdata/poisson_advanced4/sort_4/A_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_0_6_double("cuda q2 sort_0 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_0_7_double("cuda q2 sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_0_8_double("cuda q2 sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_0/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_1_6_double("cuda q2 sort_1 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_1_7_double("cuda q2 sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_1_8_double("cuda q2 sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_1/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_2_6_double("cuda q2 sort_2 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_2_7_double("cuda q2 sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_2_8_double("cuda q2 sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_2/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_3_6_double("cuda q2 sort_3 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_3_7_double("cuda q2 sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_3_8_double("cuda q2 sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_3/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_4_6_double("cuda q2 sort_4 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_4_7_double("cuda q2 sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q2_sort_4_8_double("cuda q2 sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q2_sort_4/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_0_6_double("cuda q1t sort_0 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_0_7_double("cuda q1t sort_0 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_0_8_double("cuda q1t sort_0 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_1_6_double("cuda q1t sort_1 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_1_7_double("cuda q1t sort_1 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_1_8_double("cuda q1t sort_1 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_2_6_double("cuda q1t sort_2 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_2_7_double("cuda q1t sort_2 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_2_8_double("cuda q1t sort_2 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_3_6_double("cuda q1t sort_3 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_3_7_double("cuda q1t sort_3 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_3_8_double("cuda q1t sort_3 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_4_6_double("cuda q1t sort_4 double L6 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_4_7_double("cuda q1t sort_4 double L7 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_q1t_sort_4_8_double("cuda q1t sort_4 double L8 ELLPACK stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_8.ell", 10);
//--------------------PROLMATS-----------------------------
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_0_7_double("mcsse q1 sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_0_8_double("mcsse q1 sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_0_9_double("mcsse q1 sort_0 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_1_7_double("mcsse q1 sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_1_8_double("mcsse q1 sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_1_9_double("mcsse q1 sort_1 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_2_7_double("mcsse q1 sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_2_8_double("mcsse q1 sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_2_9_double("mcsse q1 sort_2 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_3_7_double("mcsse q1 sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_3_8_double("mcsse q1 sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_3_9_double("mcsse q1 sort_3 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_4_7_double("mcsse q1 sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_4_8_double("mcsse q1 sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1_sort_4_9_double("mcsse q1 sort_4 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_0_6_double("mcsse q2 sort_0 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_0_7_double("mcsse q2 sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_0_8_double("mcsse q2 sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_1_6_double("mcsse q2 sort_1 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_1_7_double("mcsse q2 sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_1_8_double("mcsse q2 sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_2_6_double("mcsse q2 sort_2 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_2_7_double("mcsse q2 sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_2_8_double("mcsse q2 sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_3_6_double("mcsse q2 sort_3 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_3_7_double("mcsse q2 sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_3_8_double("mcsse q2 sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_4_6_double("mcsse q2 sort_4 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_4_7_double("mcsse q2 sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q2_sort_4_8_double("mcsse q2 sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_0_6_double("mcsse q1t sort_0 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_0_7_double("mcsse q1t sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_0_8_double("mcsse q1t sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_1_6_double("mcsse q1t sort_1 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_1_7_double("mcsse q1t sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_1_8_double("mcsse q1t sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_2_6_double("mcsse q1t sort_2 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_2_7_double("mcsse q1t sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_2_8_double("mcsse q1t sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_3_6_double("mcsse q1t sort_3 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_3_7_double("mcsse q1t sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_3_8_double("mcsse q1t sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_4_6_double("mcsse q1t sort_4 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_4_7_double("mcsse q1t sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_prolmat__q1t_sort_4_8_double("mcsse q1t sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_0_7_double("cuda q1 sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_0_8_double("cuda q1 sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_0_9_double("cuda q1 sort_0 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_1_7_double("cuda q1 sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_1_8_double("cuda q1 sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_1_9_double("cuda q1 sort_1 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_2_7_double("cuda q1 sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_2_8_double("cuda q1 sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_2_9_double("cuda q1 sort_2 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_3_7_double("cuda q1 sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_3_8_double("cuda q1 sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_3_9_double("cuda q1 sort_3 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_4_7_double("cuda q1 sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_4_8_double("cuda q1 sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1_sort_4_9_double("cuda q1 sort_4 double L9 ELLPACK prolmat", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_0_6_double("cuda q2 sort_0 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_0_7_double("cuda q2 sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_0_8_double("cuda q2 sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_1_6_double("cuda q2 sort_1 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_1_7_double("cuda q2 sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_1_8_double("cuda q2 sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_2_6_double("cuda q2 sort_2 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_2_7_double("cuda q2 sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_2_8_double("cuda q2 sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_3_6_double("cuda q2 sort_3 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_3_7_double("cuda q2 sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_3_8_double("cuda q2 sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_4_6_double("cuda q2 sort_4 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_4_7_double("cuda q2 sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q2_sort_4_8_double("cuda q2 sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_0_6_double("cuda q1t sort_0 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_0_7_double("cuda q1t sort_0 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_0_8_double("cuda q1t sort_0 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_1_6_double("cuda q1t sort_1 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_1_7_double("cuda q1t sort_1 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_1_8_double("cuda q1t sort_1 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_2_6_double("cuda q1t sort_2 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_2_7_double("cuda q1t sort_2 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_2_8_double("cuda q1t sort_2 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_3_6_double("cuda q1t sort_3 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_3_7_double("cuda q1t sort_3 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_3_8_double("cuda q1t sort_3 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_4_6_double("cuda q1t sort_4 double L6 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_4_7_double("cuda q1t sort_4 double L7 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_prolmat__q1t_sort_4_8_double("cuda q1t sort_4 double L8 ELLPACK prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
//--------------------Precon-----------------------------
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_0_7_double("mcsse q1 sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_7_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_0_8_double("mcsse q1 sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_8_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_0_9_double("mcsse q1 sort_0 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_9_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_1_7_double("mcsse q1 sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_7_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_1_8_double("mcsse q1 sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_8_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_1_9_double("mcsse q1 sort_1 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_9_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_2_7_double("mcsse q1 sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_7_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_2_8_double("mcsse q1 sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_8_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_2_9_double("mcsse q1 sort_2 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_9_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_3_7_double("mcsse q1 sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_7_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_3_8_double("mcsse q1 sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_8_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_3_9_double("mcsse q1 sort_3 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_9_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_4_7_double("mcsse q1 sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_7_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_4_8_double("mcsse q1 sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_8_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1_sort_4_9_double("mcsse q1 sort_4 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_9_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_0_6_double("mcsse q2 sort_0 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_0_7_double("mcsse q2 sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_0_8_double("mcsse q2 sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_1_6_double("mcsse q2 sort_1 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_1_7_double("mcsse q2 sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_1_8_double("mcsse q2 sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_2_6_double("mcsse q2 sort_2 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_2_7_double("mcsse q2 sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_2_8_double("mcsse q2 sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_3_6_double("mcsse q2 sort_3 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_3_7_double("mcsse q2 sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_3_8_double("mcsse q2 sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_4_6_double("mcsse q2 sort_4 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_4_7_double("mcsse q2 sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q2_sort_4_8_double("mcsse q2 sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_0_6_double("mcsse q1t sort_0 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_0_7_double("mcsse q1t sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_0_8_double("mcsse q1t sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_1_6_double("mcsse q1t sort_1 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_1_7_double("mcsse q1t sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_1_8_double("mcsse q1t sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_2_6_double("mcsse q1t sort_2 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_2_7_double("mcsse q1t sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_2_8_double("mcsse q1t sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_3_6_double("mcsse q1t sort_3 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_3_7_double("mcsse q1t sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_3_8_double("mcsse q1t sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_8_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_4_6_double("mcsse q1t sort_4 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_6_spai.ell", 10);
Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_4_7_double("mcsse q1t sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_7_spai.ell", 10);
//Pareng<tags::CPU::MultiCore::SSE, double> mcsse_pareng_ell_precon__q1t_sort_4_8_double("mcsse q1t sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_0_7_double("cuda q1 sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_7_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_0_8_double("cuda q1 sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_8_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_0_9_double("cuda q1 sort_0 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_0/A_9_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_1_7_double("cuda q1 sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_7_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_1_8_double("cuda q1 sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_8_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_1_9_double("cuda q1 sort_1 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_1/A_9_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_2_7_double("cuda q1 sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_7_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_2_8_double("cuda q1 sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_8_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_2_9_double("cuda q1 sort_2 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_2/A_9_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_3_7_double("cuda q1 sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_7_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_3_8_double("cuda q1 sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_8_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_3_9_double("cuda q1 sort_3 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_3/A_9_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_4_7_double("cuda q1 sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_7_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_4_8_double("cuda q1 sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_8_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1_sort_4_9_double("cuda q1 sort_4 double L9 ELLPACK precon", "testdata/poisson_advanced4/sort_4/A_9_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_0_6_double("cuda q2 sort_0 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_0_7_double("cuda q2 sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_0_8_double("cuda q2 sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_0/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_1_6_double("cuda q2 sort_1 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_1_7_double("cuda q2 sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_1_8_double("cuda q2 sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_1/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_2_6_double("cuda q2 sort_2 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_2_7_double("cuda q2 sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_2_8_double("cuda q2 sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_2/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_3_6_double("cuda q2 sort_3 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_3_7_double("cuda q2 sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_3_8_double("cuda q2 sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_3/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_4_6_double("cuda q2 sort_4 double L6 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_4_7_double("cuda q2 sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q2_sort_4_8_double("cuda q2 sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/q2_sort_4/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_0_6_double("cuda q1t sort_0 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_0_7_double("cuda q1t sort_0 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_0_8_double("cuda q1t sort_0 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_0/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_1_6_double("cuda q1t sort_1 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_1_7_double("cuda q1t sort_1 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_1_8_double("cuda q1t sort_1 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_1/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_2_6_double("cuda q1t sort_2 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_2_7_double("cuda q1t sort_2 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_2_8_double("cuda q1t sort_2 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_2/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_3_6_double("cuda q1t sort_3 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_3_7_double("cuda q1t sort_3 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_3_8_double("cuda q1t sort_3 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_3/A_8_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_4_6_double("cuda q1t sort_4 double L6 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_6_spai.ell", 10);
Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_4_7_double("cuda q1t sort_4 double L7 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_7_spai.ell", 10);
//Pareng<tags::GPU::CUDA, double> cuda_pareng_ell_precon__q1t_sort_4_8_double("cuda q1t sort_4 double L8 ELLPACK precon", "testdata/poisson_advanced4/q1t_sort_4/A_8_spai.ell", 10);
#endif
#endif

template <typename Tag_, typename DT_>
class ParengCSR:
    public Benchmark
{
    private:
        std::string _file_name;
        unsigned long _count;

    public:
        ParengCSR(const std::string & tag, std::string filename, unsigned long count) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _file_name = filename;
            _count = count;
        }

        virtual void run()
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/";
            _file_name = filebase + _file_name;
            SparseMatrixELL<DT_> smatrix = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            SparseMatrixCSR<DT_> smatrix_final(smatrix);

            DenseVector<DT_> x(smatrix_final.columns());
            DenseVector<DT_> y(smatrix_final.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 100 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix_final, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            {
            BenchmarkInfo info;
            info.flops = smatrix_final.used_elements() * 2;
            info.load = smatrix_final.used_elements() * 3 * sizeof(DT_);
            info.store = smatrix_final.used_elements() * 1 * sizeof(DT_);
            evaluate(info * 100);
            }
            std::cout<<"Non Zero Elements: "<<smatrix_final.used_elements()<<std::endl;
        }
};

#ifdef HONEI_SSE
#ifdef HONEI_CUDA_DOUBLE
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_0_7_double("mcsse q1 sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_0_8_double("mcsse q1 sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_0_9_double("mcsse q1 sort_0 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_1_7_double("mcsse q1 sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_1_8_double("mcsse q1 sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_1_9_double("mcsse q1 sort_1 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_2_7_double("mcsse q1 sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_2_8_double("mcsse q1 sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_2_9_double("mcsse q1 sort_2 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_3_7_double("mcsse q1 sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_3_8_double("mcsse q1 sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_3_9_double("mcsse q1 sort_3 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_4_7_double("mcsse q1 sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_4_8_double("mcsse q1 sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1_sort_4_9_double("mcsse q1 sort_4 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_0_6_double("mcsse q2 sort_0 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_0_7_double("mcsse q2 sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_0_8_double("mcsse q2 sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_1_6_double("mcsse q2 sort_1 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_1_7_double("mcsse q2 sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_1_8_double("mcsse q2 sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_2_6_double("mcsse q2 sort_2 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_2_7_double("mcsse q2 sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_2_8_double("mcsse q2 sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_3_6_double("mcsse q2 sort_3 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_3_7_double("mcsse q2 sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_3_8_double("mcsse q2 sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_4_6_double("mcsse q2 sort_4 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_4_7_double("mcsse q2 sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q2_sort_4_8_double("mcsse q2 sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_0_6_double("mcsse q1t sort_0 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_0_7_double("mcsse q1t sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_0_8_double("mcsse q1t sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_1_6_double("mcsse q1t sort_1 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_1_7_double("mcsse q1t sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_1_8_double("mcsse q1t sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_2_6_double("mcsse q1t sort_2 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_2_7_double("mcsse q1t sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_2_8_double("mcsse q1t sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_3_6_double("mcsse q1t sort_3 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_3_7_double("mcsse q1t sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_3_8_double("mcsse q1t sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_4_6_double("mcsse q1t sort_4 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_4_7_double("mcsse q1t sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_q1t_sort_4_8_double("mcsse q1t sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_0_7_double("cuda q1 sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_0_8_double("cuda q1 sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_0_9_double("cuda q1 sort_0 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_0/A_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_1_7_double("cuda q1 sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_1_8_double("cuda q1 sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_1_9_double("cuda q1 sort_1 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_1/A_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_2_7_double("cuda q1 sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_2_8_double("cuda q1 sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_2_9_double("cuda q1 sort_2 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_2/A_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_3_7_double("cuda q1 sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_3_8_double("cuda q1 sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_3_9_double("cuda q1 sort_3 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_3/A_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_4_7_double("cuda q1 sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_4_8_double("cuda q1 sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1_sort_4_9_double("cuda q1 sort_4 double L9 CSR stiffness", "testdata/poisson_advanced4/sort_4/A_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_0_6_double("cuda q2 sort_0 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_0_7_double("cuda q2 sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_0_8_double("cuda q2 sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_0/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_1_6_double("cuda q2 sort_1 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_1_7_double("cuda q2 sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_1_8_double("cuda q2 sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_1/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_2_6_double("cuda q2 sort_2 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_2_7_double("cuda q2 sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_2_8_double("cuda q2 sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_2/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_3_6_double("cuda q2 sort_3 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_3_7_double("cuda q2 sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_3_8_double("cuda q2 sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_3/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_4_6_double("cuda q2 sort_4 double L6 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_4_7_double("cuda q2 sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q2_sort_4_8_double("cuda q2 sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/q2_sort_4/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_0_6_double("cuda q1t sort_0 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_0_7_double("cuda q1t sort_0 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_0_8_double("cuda q1t sort_0 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_0/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_1_6_double("cuda q1t sort_1 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_1_7_double("cuda q1t sort_1 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_1_8_double("cuda q1t sort_1 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_1/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_2_6_double("cuda q1t sort_2 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_2_7_double("cuda q1t sort_2 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_2_8_double("cuda q1t sort_2 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_2/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_3_6_double("cuda q1t sort_3 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_3_7_double("cuda q1t sort_3 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_3_8_double("cuda q1t sort_3 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_3/A_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_4_6_double("cuda q1t sort_4 double L6 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_4_7_double("cuda q1t sort_4 double L7 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_q1t_sort_4_8_double("cuda q1t sort_4 double L8 CSR stiffness", "testdata/poisson_advanced4/q1t_sort_4/A_8.ell", 10);
//--------------------PROLMATS-----------------------------
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_0_7_double("mcsse q1 sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_0_8_double("mcsse q1 sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_0_9_double("mcsse q1 sort_0 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_1_7_double("mcsse q1 sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_1_8_double("mcsse q1 sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_1_9_double("mcsse q1 sort_1 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_2_7_double("mcsse q1 sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_2_8_double("mcsse q1 sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_2_9_double("mcsse q1 sort_2 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_3_7_double("mcsse q1 sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_3_8_double("mcsse q1 sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_3_9_double("mcsse q1 sort_3 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_4_7_double("mcsse q1 sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_4_8_double("mcsse q1 sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1_sort_4_9_double("mcsse q1 sort_4 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_0_6_double("mcsse q2 sort_0 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_0_7_double("mcsse q2 sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_0_8_double("mcsse q2 sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_1_6_double("mcsse q2 sort_1 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_1_7_double("mcsse q2 sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_1_8_double("mcsse q2 sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_2_6_double("mcsse q2 sort_2 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_2_7_double("mcsse q2 sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_2_8_double("mcsse q2 sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_3_6_double("mcsse q2 sort_3 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_3_7_double("mcsse q2 sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_3_8_double("mcsse q2 sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_4_6_double("mcsse q2 sort_4 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_4_7_double("mcsse q2 sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q2_sort_4_8_double("mcsse q2 sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_0_6_double("mcsse q1t sort_0 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_0_7_double("mcsse q1t sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_0_8_double("mcsse q1t sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_1_6_double("mcsse q1t sort_1 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_1_7_double("mcsse q1t sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_1_8_double("mcsse q1t sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_2_6_double("mcsse q1t sort_2 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_2_7_double("mcsse q1t sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_2_8_double("mcsse q1t sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_3_6_double("mcsse q1t sort_3 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_3_7_double("mcsse q1t sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_3_8_double("mcsse q1t sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_4_6_double("mcsse q1t sort_4 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_4_7_double("mcsse q1t sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_prolmat__q1t_sort_4_8_double("mcsse q1t sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_0_7_double("cuda q1 sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_0_8_double("cuda q1 sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_0_9_double("cuda q1 sort_0 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_1_7_double("cuda q1 sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_1_8_double("cuda q1 sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_1_9_double("cuda q1 sort_1 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_2_7_double("cuda q1 sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_2_8_double("cuda q1 sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_2_9_double("cuda q1 sort_2 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_3_7_double("cuda q1 sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_3_8_double("cuda q1 sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_3_9_double("cuda q1 sort_3 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_4_7_double("cuda q1 sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_4_8_double("cuda q1 sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1_sort_4_9_double("cuda q1 sort_4 double L9 CSR prolmat", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_0_6_double("cuda q2 sort_0 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_0_7_double("cuda q2 sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_0_8_double("cuda q2 sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_1_6_double("cuda q2 sort_1 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_1_7_double("cuda q2 sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_1_8_double("cuda q2 sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_2_6_double("cuda q2 sort_2 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_2_7_double("cuda q2 sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_2_8_double("cuda q2 sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_3_6_double("cuda q2 sort_3 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_3_7_double("cuda q2 sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_3_8_double("cuda q2 sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_4_6_double("cuda q2 sort_4 double L6 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_4_7_double("cuda q2 sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q2_sort_4_8_double("cuda q2 sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_0_6_double("cuda q1t sort_0 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_0_7_double("cuda q1t sort_0 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_0_8_double("cuda q1t sort_0 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_1_6_double("cuda q1t sort_1 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_1_7_double("cuda q1t sort_1 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_1_8_double("cuda q1t sort_1 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_2_6_double("cuda q1t sort_2 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_2_7_double("cuda q1t sort_2 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_2_8_double("cuda q1t sort_2 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_3_6_double("cuda q1t sort_3 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_3_7_double("cuda q1t sort_3 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_3_8_double("cuda q1t sort_3 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_4_6_double("cuda q1t sort_4 double L6 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_4_7_double("cuda q1t sort_4 double L7 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_prolmat__q1t_sort_4_8_double("cuda q1t sort_4 double L8 CSR prolmat", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
//--------------------Precon-----------------------------
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_0_7_double("mcsse q1 sort_0 double L7 CSR precon", "testdata/poisson_advanced4/sort_0/A_7_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_0_8_double("mcsse q1 sort_0 double L8 CSR precon", "testdata/poisson_advanced4/sort_0/A_8_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_0_9_double("mcsse q1 sort_0 double L9 CSR precon", "testdata/poisson_advanced4/sort_0/A_9_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_1_7_double("mcsse q1 sort_1 double L7 CSR precon", "testdata/poisson_advanced4/sort_1/A_7_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_1_8_double("mcsse q1 sort_1 double L8 CSR precon", "testdata/poisson_advanced4/sort_1/A_8_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_1_9_double("mcsse q1 sort_1 double L9 CSR precon", "testdata/poisson_advanced4/sort_1/A_9_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_2_7_double("mcsse q1 sort_2 double L7 CSR precon", "testdata/poisson_advanced4/sort_2/A_7_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_2_8_double("mcsse q1 sort_2 double L8 CSR precon", "testdata/poisson_advanced4/sort_2/A_8_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_2_9_double("mcsse q1 sort_2 double L9 CSR precon", "testdata/poisson_advanced4/sort_2/A_9_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_3_7_double("mcsse q1 sort_3 double L7 CSR precon", "testdata/poisson_advanced4/sort_3/A_7_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_3_8_double("mcsse q1 sort_3 double L8 CSR precon", "testdata/poisson_advanced4/sort_3/A_8_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_3_9_double("mcsse q1 sort_3 double L9 CSR precon", "testdata/poisson_advanced4/sort_3/A_9_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_4_7_double("mcsse q1 sort_4 double L7 CSR precon", "testdata/poisson_advanced4/sort_4/A_7_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_4_8_double("mcsse q1 sort_4 double L8 CSR precon", "testdata/poisson_advanced4/sort_4/A_8_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1_sort_4_9_double("mcsse q1 sort_4 double L9 CSR precon", "testdata/poisson_advanced4/sort_4/A_9_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_0_6_double("mcsse q2 sort_0 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_0_7_double("mcsse q2 sort_0 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_0_8_double("mcsse q2 sort_0 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_1_6_double("mcsse q2 sort_1 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_1_7_double("mcsse q2 sort_1 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_1_8_double("mcsse q2 sort_1 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_2_6_double("mcsse q2 sort_2 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_2_7_double("mcsse q2 sort_2 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_2_8_double("mcsse q2 sort_2 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_3_6_double("mcsse q2 sort_3 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_3_7_double("mcsse q2 sort_3 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_3_8_double("mcsse q2 sort_3 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_4_6_double("mcsse q2 sort_4 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_4_7_double("mcsse q2 sort_4 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q2_sort_4_8_double("mcsse q2 sort_4 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_0_6_double("mcsse q1t sort_0 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_0_7_double("mcsse q1t sort_0 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_0_8_double("mcsse q1t sort_0 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_1_6_double("mcsse q1t sort_1 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_1_7_double("mcsse q1t sort_1 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_1_8_double("mcsse q1t sort_1 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_2_6_double("mcsse q1t sort_2 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_2_7_double("mcsse q1t sort_2 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_2_8_double("mcsse q1t sort_2 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_3_6_double("mcsse q1t sort_3 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_3_7_double("mcsse q1t sort_3 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_3_8_double("mcsse q1t sort_3 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_8_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_4_6_double("mcsse q1t sort_4 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_6_spai.ell", 10);
ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_4_7_double("mcsse q1t sort_4 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_7_spai.ell", 10);
//ParengCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_CSR_precon__q1t_sort_4_8_double("mcsse q1t sort_4 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_0_7_double("cuda q1 sort_0 double L7 CSR precon", "testdata/poisson_advanced4/sort_0/A_7_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_0_8_double("cuda q1 sort_0 double L8 CSR precon", "testdata/poisson_advanced4/sort_0/A_8_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_0_9_double("cuda q1 sort_0 double L9 CSR precon", "testdata/poisson_advanced4/sort_0/A_9_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_1_7_double("cuda q1 sort_1 double L7 CSR precon", "testdata/poisson_advanced4/sort_1/A_7_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_1_8_double("cuda q1 sort_1 double L8 CSR precon", "testdata/poisson_advanced4/sort_1/A_8_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_1_9_double("cuda q1 sort_1 double L9 CSR precon", "testdata/poisson_advanced4/sort_1/A_9_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_2_7_double("cuda q1 sort_2 double L7 CSR precon", "testdata/poisson_advanced4/sort_2/A_7_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_2_8_double("cuda q1 sort_2 double L8 CSR precon", "testdata/poisson_advanced4/sort_2/A_8_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_2_9_double("cuda q1 sort_2 double L9 CSR precon", "testdata/poisson_advanced4/sort_2/A_9_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_3_7_double("cuda q1 sort_3 double L7 CSR precon", "testdata/poisson_advanced4/sort_3/A_7_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_3_8_double("cuda q1 sort_3 double L8 CSR precon", "testdata/poisson_advanced4/sort_3/A_8_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_3_9_double("cuda q1 sort_3 double L9 CSR precon", "testdata/poisson_advanced4/sort_3/A_9_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_4_7_double("cuda q1 sort_4 double L7 CSR precon", "testdata/poisson_advanced4/sort_4/A_7_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_4_8_double("cuda q1 sort_4 double L8 CSR precon", "testdata/poisson_advanced4/sort_4/A_8_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1_sort_4_9_double("cuda q1 sort_4 double L9 CSR precon", "testdata/poisson_advanced4/sort_4/A_9_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_0_6_double("cuda q2 sort_0 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_0_7_double("cuda q2 sort_0 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_0_8_double("cuda q2 sort_0 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_0/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_1_6_double("cuda q2 sort_1 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_1_7_double("cuda q2 sort_1 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_1_8_double("cuda q2 sort_1 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_1/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_2_6_double("cuda q2 sort_2 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_2_7_double("cuda q2 sort_2 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_2_8_double("cuda q2 sort_2 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_2/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_3_6_double("cuda q2 sort_3 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_3_7_double("cuda q2 sort_3 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_3_8_double("cuda q2 sort_3 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_3/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_4_6_double("cuda q2 sort_4 double L6 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_4_7_double("cuda q2 sort_4 double L7 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q2_sort_4_8_double("cuda q2 sort_4 double L8 CSR precon", "testdata/poisson_advanced4/q2_sort_4/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_0_6_double("cuda q1t sort_0 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_0_7_double("cuda q1t sort_0 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_0_8_double("cuda q1t sort_0 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_0/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_1_6_double("cuda q1t sort_1 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_1_7_double("cuda q1t sort_1 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_1_8_double("cuda q1t sort_1 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_1/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_2_6_double("cuda q1t sort_2 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_2_7_double("cuda q1t sort_2 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_2_8_double("cuda q1t sort_2 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_2/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_3_6_double("cuda q1t sort_3 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_3_7_double("cuda q1t sort_3 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_3_8_double("cuda q1t sort_3 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_3/A_8_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_4_6_double("cuda q1t sort_4 double L6 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_6_spai.ell", 10);
ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_4_7_double("cuda q1t sort_4 double L7 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_7_spai.ell", 10);
//ParengCSR<tags::GPU::CUDA, double> cuda_pareng_CSR_precon__q1t_sort_4_8_double("cuda q1t sort_4 double L8 CSR precon", "testdata/poisson_advanced4/q1t_sort_4/A_8_spai.ell", 10);
#endif
#endif

template <typename Tag_, typename DT_>
class ParengResmat:
    public Benchmark
{
    private:
        std::string _file_name;
        unsigned long _count;

    public:
        ParengResmat(const std::string & tag, std::string filename, unsigned long count) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _file_name = filename;
            _count = count;
        }

        virtual void run()
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/";
            _file_name = filebase + _file_name;

            SparseMatrixELL<DT_> smatrix_temp = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            SparseMatrix<DT_> smatrix_temp_sparse(smatrix_temp);
            SparseMatrix<DT_> smatrix_temp2(smatrix_temp.columns(), smatrix_temp.rows());
            Transposition<tags::CPU>::value(smatrix_temp_sparse, smatrix_temp2);
            SparseMatrixELL<DT_> smatrix(smatrix_temp2);

            DenseVector<DT_> x(smatrix.columns());
            DenseVector<DT_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 100 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            {
            BenchmarkInfo info;
            info.flops = smatrix.used_elements() * 2;
            info.load = smatrix.used_elements() * 3 * sizeof(DT_);
            info.store = smatrix.used_elements() * 1 * sizeof(DT_);
            evaluate(info * 100);
            }
            std::cout<<"Non Zero Elements: "<<smatrix.used_elements()<<std::endl;
        }
};
#ifdef HONEI_SSE
#ifdef HONEI_CUDA_DOUBLE
//--------------------Resmat ELLPACKELL-----------------------------
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_0_7_double("mcsse q1 sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_0_8_double("mcsse q1 sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_0_9_double("mcsse q1 sort_0 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_1_7_double("mcsse q1 sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_1_8_double("mcsse q1 sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_1_9_double("mcsse q1 sort_1 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_2_7_double("mcsse q1 sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_2_8_double("mcsse q1 sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_2_9_double("mcsse q1 sort_2 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_3_7_double("mcsse q1 sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_3_8_double("mcsse q1 sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_3_9_double("mcsse q1 sort_3 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_4_7_double("mcsse q1 sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_4_8_double("mcsse q1 sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1_sort_4_9_double("mcsse q1 sort_4 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_0_6_double("mcsse q2 sort_0 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_0_7_double("mcsse q2 sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_0_8_double("mcsse q2 sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_1_6_double("mcsse q2 sort_1 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_1_7_double("mcsse q2 sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_1_8_double("mcsse q2 sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_2_6_double("mcsse q2 sort_2 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_2_7_double("mcsse q2 sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_2_8_double("mcsse q2 sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_3_6_double("mcsse q2 sort_3 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_3_7_double("mcsse q2 sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_3_8_double("mcsse q2 sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_4_6_double("mcsse q2 sort_4 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_4_7_double("mcsse q2 sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q2_sort_4_8_double("mcsse q2 sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_0_6_double("mcsse q1t sort_0 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_0_7_double("mcsse q1t sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_0_8_double("mcsse q1t sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_1_6_double("mcsse q1t sort_1 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_1_7_double("mcsse q1t sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_1_8_double("mcsse q1t sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_2_6_double("mcsse q1t sort_2 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_2_7_double("mcsse q1t sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_2_8_double("mcsse q1t sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_3_6_double("mcsse q1t sort_3 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_3_7_double("mcsse q1t sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_3_8_double("mcsse q1t sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_4_6_double("mcsse q1t sort_4 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_4_7_double("mcsse q1t sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengResmat<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat___q1t_sort_4_8_double("mcsse q1t sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_0_7_double("cuda q1 sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_0_8_double("cuda q1 sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_0_9_double("cuda q1 sort_0 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_1_7_double("cuda q1 sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_1_8_double("cuda q1 sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_1_9_double("cuda q1 sort_1 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_2_7_double("cuda q1 sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_2_8_double("cuda q1 sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_2_9_double("cuda q1 sort_2 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_3_7_double("cuda q1 sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_3_8_double("cuda q1 sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_3_9_double("cuda q1 sort_3 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_4_7_double("cuda q1 sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_4_8_double("cuda q1 sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1_sort_4_9_double("cuda q1 sort_4 double L9 Resmat ELLPACK", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_0_6_double("cuda q2 sort_0 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_0_7_double("cuda q2 sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_0_8_double("cuda q2 sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_1_6_double("cuda q2 sort_1 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_1_7_double("cuda q2 sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_1_8_double("cuda q2 sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_2_6_double("cuda q2 sort_2 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_2_7_double("cuda q2 sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_2_8_double("cuda q2 sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_3_6_double("cuda q2 sort_3 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_3_7_double("cuda q2 sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_3_8_double("cuda q2 sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_4_6_double("cuda q2 sort_4 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_4_7_double("cuda q2 sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q2_sort_4_8_double("cuda q2 sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_0_6_double("cuda q1t sort_0 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_0_7_double("cuda q1t sort_0 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_0_8_double("cuda q1t sort_0 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_1_6_double("cuda q1t sort_1 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_1_7_double("cuda q1t sort_1 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_1_8_double("cuda q1t sort_1 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_2_6_double("cuda q1t sort_2 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_2_7_double("cuda q1t sort_2 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_2_8_double("cuda q1t sort_2 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_3_6_double("cuda q1t sort_3 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_3_7_double("cuda q1t sort_3 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_3_8_double("cuda q1t sort_3 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_4_6_double("cuda q1t sort_4 double L6 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_4_7_double("cuda q1t sort_4 double L7 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengResmat<tags::GPU::CUDA, double> cuda_pareng_Resmat___q1t_sort_4_8_double("cuda q1t sort_4 double L8 Resmat ELLPACK", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
#endif
#endif

template <typename Tag_, typename DT_>
class ParengResmatCSR:
    public Benchmark
{
    private:
        std::string _file_name;
        unsigned long _count;

    public:
        ParengResmatCSR(const std::string & tag, std::string filename, unsigned long count) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _file_name = filename;
            _count = count;
        }

        virtual void run()
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/";
            _file_name = filebase + _file_name;

            SparseMatrixELL<DT_> smatrix_temp = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            SparseMatrix<DT_> smatrix_temp_sparse(smatrix_temp);
            SparseMatrix<DT_> smatrix_temp2(smatrix_temp.columns(), smatrix_temp.rows());
            Transposition<tags::CPU>::value(smatrix_temp_sparse, smatrix_temp2);
            SparseMatrixCSR<DT_> smatrix(smatrix_temp2);

            DenseVector<DT_> x(smatrix.columns());
            DenseVector<DT_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 100 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            {
            BenchmarkInfo info;
            info.flops = smatrix.used_elements() * 2;
            info.load = smatrix.used_elements() * 3 * sizeof(DT_);
            info.store = smatrix.used_elements() * 1 * sizeof(DT_);
            evaluate(info * 100);
            }
            std::cout<<"Non Zero Elements: "<<smatrix.used_elements()<<std::endl;
        }
};
#ifdef HONEI_SSE
#ifdef HONEI_CUDA_DOUBLE
//--------------------Resmat CSRELL-----------------------------
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_0_7_double("mcsse q1 sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_0_8_double("mcsse q1 sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_0_9_double("mcsse q1 sort_0 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_1_7_double("mcsse q1 sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_1_8_double("mcsse q1 sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_1_9_double("mcsse q1 sort_1 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_2_7_double("mcsse q1 sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_2_8_double("mcsse q1 sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_2_9_double("mcsse q1 sort_2 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_3_7_double("mcsse q1 sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_3_8_double("mcsse q1 sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_3_9_double("mcsse q1 sort_3 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_4_7_double("mcsse q1 sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_4_8_double("mcsse q1 sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1_sort_4_9_double("mcsse q1 sort_4 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_0_6_double("mcsse q2 sort_0 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_0_7_double("mcsse q2 sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_0_8_double("mcsse q2 sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_1_6_double("mcsse q2 sort_1 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_1_7_double("mcsse q2 sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_1_8_double("mcsse q2 sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_2_6_double("mcsse q2 sort_2 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_2_7_double("mcsse q2 sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_2_8_double("mcsse q2 sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_3_6_double("mcsse q2 sort_3 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_3_7_double("mcsse q2 sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_3_8_double("mcsse q2 sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_4_6_double("mcsse q2 sort_4 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_4_7_double("mcsse q2 sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q2_sort_4_8_double("mcsse q2 sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_0_6_double("mcsse q1t sort_0 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_0_7_double("mcsse q1t sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_0_8_double("mcsse q1t sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_1_6_double("mcsse q1t sort_1 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_1_7_double("mcsse q1t sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_1_8_double("mcsse q1t sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_2_6_double("mcsse q1t sort_2 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_2_7_double("mcsse q1t sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_2_8_double("mcsse q1t sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_3_6_double("mcsse q1t sort_3 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_3_7_double("mcsse q1t sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_3_8_double("mcsse q1t sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_4_6_double("mcsse q1t sort_4 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_4_7_double("mcsse q1t sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengResmatCSR<tags::CPU::MultiCore::SSE, double> mcsse_pareng_Resmat_CSR_q1t_sort_4_8_double("mcsse q1t sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_0_7_double("cuda q1 sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_7.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_0_8_double("cuda q1 sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_8.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_0_9_double("cuda q1 sort_0 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_0/prol_9.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_1_7_double("cuda q1 sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_7.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_1_8_double("cuda q1 sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_8.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_1_9_double("cuda q1 sort_1 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_1/prol_9.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_2_7_double("cuda q1 sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_7.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_2_8_double("cuda q1 sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_8.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_2_9_double("cuda q1 sort_2 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_2/prol_9.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_3_7_double("cuda q1 sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_7.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_3_8_double("cuda q1 sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_8.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_3_9_double("cuda q1 sort_3 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_3/prol_9.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_4_7_double("cuda q1 sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_7.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_4_8_double("cuda q1 sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_8.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1_sort_4_9_double("cuda q1 sort_4 double L9 Resmat CSR", "testdata/poisson_advanced4/sort_4/prol_9.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_0_6_double("cuda q2 sort_0 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_0_7_double("cuda q2 sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_0_8_double("cuda q2 sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_0/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_1_6_double("cuda q2 sort_1 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_1_7_double("cuda q2 sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_1_8_double("cuda q2 sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_1/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_2_6_double("cuda q2 sort_2 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_2_7_double("cuda q2 sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_2_8_double("cuda q2 sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_2/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_3_6_double("cuda q2 sort_3 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_3_7_double("cuda q2 sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_3_8_double("cuda q2 sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_3/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_4_6_double("cuda q2 sort_4 double L6 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_4_7_double("cuda q2 sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q2_sort_4_8_double("cuda q2 sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/q2_sort_4/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_0_6_double("cuda q1t sort_0 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_0_7_double("cuda q1t sort_0 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_0_8_double("cuda q1t sort_0 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_0/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_1_6_double("cuda q1t sort_1 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_1_7_double("cuda q1t sort_1 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_1_8_double("cuda q1t sort_1 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_1/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_2_6_double("cuda q1t sort_2 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_2_7_double("cuda q1t sort_2 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_2_8_double("cuda q1t sort_2 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_2/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_3_6_double("cuda q1t sort_3 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_3_7_double("cuda q1t sort_3 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_3_8_double("cuda q1t sort_3 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_3/prol_8.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_4_6_double("cuda q1t sort_4 double L6 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_6.ell", 10);
ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_4_7_double("cuda q1t sort_4 double L7 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_7.ell", 10);
//ParengResmatCSR<tags::GPU::CUDA, double> cuda_pareng_Resmat_CSR_q1t_sort_4_8_double("cuda q1t sort_4 double L8 Resmat CSR", "testdata/poisson_advanced4/q1t_sort_4/prol_8.ell", 10);
#endif
#endif
