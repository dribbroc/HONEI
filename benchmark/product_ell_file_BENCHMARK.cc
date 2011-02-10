/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmmund.de>
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
#include <benchmark/benchmark.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/backends/cuda/gpu_pool.hh>
#ifdef HONEI_OPENCL
#include <honei/backends/opencl/opencl_backend.hh>
#endif

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace std;


template <typename Tag_, typename DT_>
class ProductELLFileBenchmark:
    public Benchmark
{
    private:
        std::string _file_name;
        unsigned long _count;

    public:
        ProductELLFileBenchmark(const std::string & tag, std::string filename, unsigned long count) :
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
            SparseMatrixELL<DT_> bsmatrix = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            SparseMatrix<DT_> bla(bsmatrix);
            SparseMatrixELL<DT_>smatrix(bla, 1);
            std::cout<<smatrix.num_cols_per_row()<<" "<<smatrix.rows()<<std::endl;
            DenseVector<DT_> x(smatrix.rows());
            DenseVector<DT_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
#ifdef HONEI_OPENCL
                        if (Tag_::tag_value == tags::tv_opencl)
                            OpenCLBackend::instance()->flush();
#endif
                        );
            }
            {
            BenchmarkInfo info;
            info.flops = smatrix.used_elements() * 2;
            info.load = smatrix.used_elements() * 3 * sizeof(DT_);
            info.store = smatrix.used_elements() * 1 * sizeof(DT_);
            evaluate(info * 10);
            }
            std::cout<<"Non Zero Elements: "<<smatrix.used_elements()<<std::endl;
        }
};
#ifdef HONEI_SSE
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q1_0("ELL  Product double sse L7, q1 sort 0", "testdata/poisson_advanced/sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q1_0("ELL  Product double sse L8, q1 sort 0", "testdata/poisson_advanced/sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q1_0("ELL  Product double sse L9, q1 sort 0", "testdata/poisson_advanced/sort_0/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_10_double_q1_0("ELL  Product double sse L10, q1 sort 0", "testdata/poisson_advanced/sort_0/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q1_1("ELL  Product double sse L7, q1 sort 1", "testdata/poisson_advanced/sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q1_1("ELL  Product double sse L8, q1 sort 1", "testdata/poisson_advanced/sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q1_1("ELL  Product double sse L9, q1 sort 1", "testdata/poisson_advanced/sort_1/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_10_double_q1_1("ELL  Product double sse L10, q1 sort 1", "testdata/poisson_advanced/sort_1/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q1_2("ELL  Product double sse L7, q1 sort 2", "testdata/poisson_advanced/sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q1_2("ELL  Product double sse L8, q1 sort 2", "testdata/poisson_advanced/sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q1_2("ELL  Product double sse L9, q1 sort 2", "testdata/poisson_advanced/sort_2/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_10_double_q1_2("ELL  Product double sse L10, q1 sort 2", "testdata/poisson_advanced/sort_2/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q1_3("ELL  Product double sse L7, q1 sort 3", "testdata/poisson_advanced/sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q1_3("ELL  Product double sse L8, q1 sort 3", "testdata/poisson_advanced/sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q1_3("ELL  Product double sse L9, q1 sort 3", "testdata/poisson_advanced/sort_3/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_10_double_q1_3("ELL  Product double sse L10, q1 sort 3", "testdata/poisson_advanced/sort_3/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q1_4("ELL  Product double sse L7, q1 sort 4", "testdata/poisson_advanced/sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q1_4("ELL  Product double sse L8, q1 sort 4", "testdata/poisson_advanced/sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q1_4("ELL  Product double sse L9, q1 sort 4", "testdata/poisson_advanced/sort_4/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_10_double_q1_4("ELL  Product double sse L10, q1 sort 4", "testdata/poisson_advanced/sort_4/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_6_double_q2_0("ELL  Product double sse L6, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q2_0("ELL  Product double sse L7, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q2_0("ELL  Product double sse L8, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q2_0("ELL  Product double sse L9, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_6_double_q2_1("ELL  Product double sse L6, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q2_1("ELL  Product double sse L7, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q2_1("ELL  Product double sse L8, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q2_1("ELL  Product double sse L9, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_6_double_q2_2("ELL  Product double sse L6, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q2_2("ELL  Product double sse L7, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q2_2("ELL  Product double sse L8, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q2_2("ELL  Product double sse L9, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_6_double_q2_3("ELL  Product double sse L6, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q2_3("ELL  Product double sse L7, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q2_3("ELL  Product double sse L8, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q2_3("ELL  Product double sse L9, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_6_double_q2_4("ELL  Product double sse L6, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_7_double_q2_4("ELL  Product double sse L7, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_8_double_q2_4("ELL  Product double sse L8, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::SSE, double> sse_pareng_9_double_q2_4("ELL  Product double sse L9, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_9.ell", 10);
//------------------------

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q1_0("ELL  Product double mcsse L7, q1 sort 0", "testdata/poisson_advanced/sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q1_0("ELL  Product double mcsse L8, q1 sort 0", "testdata/poisson_advanced/sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q1_0("ELL  Product double mcsse L9, q1 sort 0", "testdata/poisson_advanced/sort_0/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_10_double_q1_0("ELL  Product double mcsse L10, q1 sort 0", "testdata/poisson_advanced/sort_0/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q1_1("ELL  Product double mcsse L7, q1 sort 1", "testdata/poisson_advanced/sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q1_1("ELL  Product double mcsse L8, q1 sort 1", "testdata/poisson_advanced/sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q1_1("ELL  Product double mcsse L9, q1 sort 1", "testdata/poisson_advanced/sort_1/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_10_double_q1_1("ELL  Product double mcsse L10, q1 sort 1", "testdata/poisson_advanced/sort_1/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q1_2("ELL  Product double mcsse L7, q1 sort 2", "testdata/poisson_advanced/sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q1_2("ELL  Product double mcsse L8, q1 sort 2", "testdata/poisson_advanced/sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q1_2("ELL  Product double mcsse L9, q1 sort 2", "testdata/poisson_advanced/sort_2/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_10_double_q1_2("ELL  Product double mcsse L10, q1 sort 2", "testdata/poisson_advanced/sort_2/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q1_3("ELL  Product double mcsse L7, q1 sort 3", "testdata/poisson_advanced/sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q1_3("ELL  Product double mcsse L8, q1 sort 3", "testdata/poisson_advanced/sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q1_3("ELL  Product double mcsse L9, q1 sort 3", "testdata/poisson_advanced/sort_3/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_10_double_q1_3("ELL  Product double mcsse L10, q1 sort 3", "testdata/poisson_advanced/sort_3/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q1_4("ELL  Product double mcsse L7, q1 sort 4", "testdata/poisson_advanced/sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q1_4("ELL  Product double mcsse L8, q1 sort 4", "testdata/poisson_advanced/sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q1_4("ELL  Product double mcsse L9, q1 sort 4", "testdata/poisson_advanced/sort_4/A_9.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_10_double_q1_4("ELL  Product double mcsse L10, q1 sort 4", "testdata/poisson_advanced/sort_4/A_10.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_6_double_q2_0("ELL  Product double mcsse L6, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q2_0("ELL  Product double mcsse L7, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q2_0("ELL  Product double mcsse L8, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q2_0("ELL  Product double mcsse L9, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_6_double_q2_1("ELL  Product double mcsse L6, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q2_1("ELL  Product double mcsse L7, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q2_1("ELL  Product double mcsse L8, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q2_1("ELL  Product double mcsse L9, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_6_double_q2_2("ELL  Product double mcsse L6, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q2_2("ELL  Product double mcsse L7, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q2_2("ELL  Product double mcsse L8, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q2_2("ELL  Product double mcsse L9, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_6_double_q2_3("ELL  Product double mcsse L6, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q2_3("ELL  Product double mcsse L7, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q2_3("ELL  Product double mcsse L8, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q2_3("ELL  Product double mcsse L9, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_9.ell", 10);

ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_6_double_q2_4("ELL  Product double mcsse L6, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_6.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_7_double_q2_4("ELL  Product double mcsse L7, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_8_double_q2_4("ELL  Product double mcsse L8, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::CPU::MultiCore::SSE, double> mcsse_pareng_9_double_q2_4("ELL  Product double mcsse L9, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_9.ell", 10);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q1_0("ELL  Product double cuda L7, q1 sort 0", "testdata/poisson_advanced/sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q1_0("ELL  Product double cuda L8, q1 sort 0", "testdata/poisson_advanced/sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q1_0("ELL  Product double cuda L9, q1 sort 0", "testdata/poisson_advanced/sort_0/A_9.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_10_double_q1_0("ELL  Product double cuda L10, q1 sort 0", "testdata/poisson_advanced/sort_0/A_10.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q1_1("ELL  Product double cuda L7, q1 sort 1", "testdata/poisson_advanced/sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q1_1("ELL  Product double cuda L8, q1 sort 1", "testdata/poisson_advanced/sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q1_1("ELL  Product double cuda L9, q1 sort 1", "testdata/poisson_advanced/sort_1/A_9.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_10_double_q1_1("ELL  Product double cuda L10, q1 sort 1", "testdata/poisson_advanced/sort_1/A_10.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q1_2("ELL  Product double cuda L7, q1 sort 2", "testdata/poisson_advanced/sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q1_2("ELL  Product double cuda L8, q1 sort 2", "testdata/poisson_advanced/sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q1_2("ELL  Product double cuda L9, q1 sort 2", "testdata/poisson_advanced/sort_2/A_9.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_10_double_q1_2("ELL  Product double cuda L10, q1 sort 2", "testdata/poisson_advanced/sort_2/A_10.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q1_3("ELL  Product double cuda L7, q1 sort 3", "testdata/poisson_advanced/sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q1_3("ELL  Product double cuda L8, q1 sort 3", "testdata/poisson_advanced/sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q1_3("ELL  Product double cuda L9, q1 sort 3", "testdata/poisson_advanced/sort_3/A_9.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_10_double_q1_3("ELL  Product double cuda L10, q1 sort 3", "testdata/poisson_advanced/sort_3/A_10.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q1_4("ELL  Product double cuda L7, q1 sort 4", "testdata/poisson_advanced/sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q1_4("ELL  Product double cuda L8, q1 sort 4", "testdata/poisson_advanced/sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q1_4("ELL  Product double cuda L9, q1 sort 4", "testdata/poisson_advanced/sort_4/A_9.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_10_double_q1_4("ELL  Product double cuda L10, q1 sort 4", "testdata/poisson_advanced/sort_4/A_10.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_6_double_q2_0("ELL  Product double cuda L6, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_6.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q2_0("ELL  Product double cuda L7, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q2_0("ELL  Product double cuda L8, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q2_0("ELL  Product double cuda L9, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_9.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_6_double_q2_1("ELL  Product double cuda L6, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_6.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q2_1("ELL  Product double cuda L7, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q2_1("ELL  Product double cuda L8, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q2_1("ELL  Product double cuda L9, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_9.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_6_double_q2_2("ELL  Product double cuda L6, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_6.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q2_2("ELL  Product double cuda L7, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q2_2("ELL  Product double cuda L8, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q2_2("ELL  Product double cuda L9, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_9.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_6_double_q2_3("ELL  Product double cuda L6, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_6.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q2_3("ELL  Product double cuda L7, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q2_3("ELL  Product double cuda L8, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q2_3("ELL  Product double cuda L9, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_9.ell", 10);

ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_6_double_q2_4("ELL  Product double cuda L6, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_6.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_7_double_q2_4("ELL  Product double cuda L7, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_8_double_q2_4("ELL  Product double cuda L8, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::GPU::CUDA, double> cuda_pareng_9_double_q2_4("ELL  Product double cuda L9, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_9.ell", 10);
#endif
#endif

#ifdef HONEI_OPENCL
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q1_0("ELL  Product double opencl_cpu L7, q1 sort 0", "testdata/poisson_advanced/sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q1_0("ELL  Product double opencl_cpu L8, q1 sort 0", "testdata/poisson_advanced/sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q1_0("ELL  Product double opencl_cpu L9, q1 sort 0", "testdata/poisson_advanced/sort_0/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_10_double_q1_0("ELL  Product double opencl_cpu L10, q1 sort 0", "testdata/poisson_advanced/sort_0/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q1_1("ELL  Product double opencl_cpu L7, q1 sort 1", "testdata/poisson_advanced/sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q1_1("ELL  Product double opencl_cpu L8, q1 sort 1", "testdata/poisson_advanced/sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q1_1("ELL  Product double opencl_cpu L9, q1 sort 1", "testdata/poisson_advanced/sort_1/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_10_double_q1_1("ELL  Product double opencl_cpu L10, q1 sort 1", "testdata/poisson_advanced/sort_1/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q1_2("ELL  Product double opencl_cpu L7, q1 sort 2", "testdata/poisson_advanced/sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q1_2("ELL  Product double opencl_cpu L8, q1 sort 2", "testdata/poisson_advanced/sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q1_2("ELL  Product double opencl_cpu L9, q1 sort 2", "testdata/poisson_advanced/sort_2/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_10_double_q1_2("ELL  Product double opencl_cpu L10, q1 sort 2", "testdata/poisson_advanced/sort_2/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q1_3("ELL  Product double opencl_cpu L7, q1 sort 3", "testdata/poisson_advanced/sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q1_3("ELL  Product double opencl_cpu L8, q1 sort 3", "testdata/poisson_advanced/sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q1_3("ELL  Product double opencl_cpu L9, q1 sort 3", "testdata/poisson_advanced/sort_3/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_10_double_q1_3("ELL  Product double opencl_cpu L10, q1 sort 3", "testdata/poisson_advanced/sort_3/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q1_4("ELL  Product double opencl_cpu L7, q1 sort 4", "testdata/poisson_advanced/sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q1_4("ELL  Product double opencl_cpu L8, q1 sort 4", "testdata/poisson_advanced/sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q1_4("ELL  Product double opencl_cpu L9, q1 sort 4", "testdata/poisson_advanced/sort_4/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_10_double_q1_4("ELL  Product double opencl_cpu L10, q1 sort 4", "testdata/poisson_advanced/sort_4/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, float> opencl_cpu_pareng_6_float_q2_0("ELL  Product float opencl_cpu L6, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q2_0("ELL  Product double opencl_cpu L7, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q2_0("ELL  Product double opencl_cpu L8, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q2_0("ELL  Product double opencl_cpu L9, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_6_double_q2_1("ELL  Product double opencl_cpu L6, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q2_1("ELL  Product double opencl_cpu L7, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q2_1("ELL  Product double opencl_cpu L8, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q2_1("ELL  Product double opencl_cpu L9, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_6_double_q2_2("ELL  Product double opencl_cpu L6, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q2_2("ELL  Product double opencl_cpu L7, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q2_2("ELL  Product double opencl_cpu L8, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q2_2("ELL  Product double opencl_cpu L9, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_6_double_q2_3("ELL  Product double opencl_cpu L6, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q2_3("ELL  Product double opencl_cpu L7, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q2_3("ELL  Product double opencl_cpu L8, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q2_3("ELL  Product double opencl_cpu L9, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_6_double_q2_4("ELL  Product double opencl_cpu L6, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_7_double_q2_4("ELL  Product double opencl_cpu L7, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_8_double_q2_4("ELL  Product double opencl_cpu L8, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::CPU, double> opencl_cpu_pareng_9_double_q2_4("ELL  Product double opencl_cpu L9, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_9.ell", 10);
#ifdef HONEI_CUDA_DOUBLE
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q1_0("ELL  Product double opencl_gpu L7, q1 sort 0", "testdata/poisson_advanced/sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q1_0("ELL  Product double opencl_gpu L8, q1 sort 0", "testdata/poisson_advanced/sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q1_0("ELL  Product double opencl_gpu L9, q1 sort 0", "testdata/poisson_advanced/sort_0/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_10_double_q1_0("ELL  Product double opencl_gpu L10, q1 sort 0", "testdata/poisson_advanced/sort_0/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q1_1("ELL  Product double opencl_gpu L7, q1 sort 1", "testdata/poisson_advanced/sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q1_1("ELL  Product double opencl_gpu L8, q1 sort 1", "testdata/poisson_advanced/sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q1_1("ELL  Product double opencl_gpu L9, q1 sort 1", "testdata/poisson_advanced/sort_1/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_10_double_q1_1("ELL  Product double opencl_gpu L10, q1 sort 1", "testdata/poisson_advanced/sort_1/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q1_2("ELL  Product double opencl_gpu L7, q1 sort 2", "testdata/poisson_advanced/sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q1_2("ELL  Product double opencl_gpu L8, q1 sort 2", "testdata/poisson_advanced/sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q1_2("ELL  Product double opencl_gpu L9, q1 sort 2", "testdata/poisson_advanced/sort_2/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_10_double_q1_2("ELL  Product double opencl_gpu L10, q1 sort 2", "testdata/poisson_advanced/sort_2/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q1_3("ELL  Product double opencl_gpu L7, q1 sort 3", "testdata/poisson_advanced/sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q1_3("ELL  Product double opencl_gpu L8, q1 sort 3", "testdata/poisson_advanced/sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q1_3("ELL  Product double opencl_gpu L9, q1 sort 3", "testdata/poisson_advanced/sort_3/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_10_double_q1_3("ELL  Product double opencl_gpu L10, q1 sort 3", "testdata/poisson_advanced/sort_3/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q1_4("ELL  Product double opencl_gpu L7, q1 sort 4", "testdata/poisson_advanced/sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q1_4("ELL  Product double opencl_gpu L8, q1 sort 4", "testdata/poisson_advanced/sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q1_4("ELL  Product double opencl_gpu L9, q1 sort 4", "testdata/poisson_advanced/sort_4/A_9.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_10_double_q1_4("ELL  Product double opencl_gpu L10, q1 sort 4", "testdata/poisson_advanced/sort_4/A_10.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, float> opencl_gpu_pareng_6_float_q2_0("ELL  Product float opencl_gpu L6, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q2_0("ELL  Product double opencl_gpu L7, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q2_0("ELL  Product double opencl_gpu L8, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q2_0("ELL  Product double opencl_gpu L9, q2 sort 0", "testdata/poisson_advanced/q2_sort_0/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_6_double_q2_1("ELL  Product double opencl_gpu L6, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q2_1("ELL  Product double opencl_gpu L7, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q2_1("ELL  Product double opencl_gpu L8, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q2_1("ELL  Product double opencl_gpu L9, q2 sort 1", "testdata/poisson_advanced/q2_sort_1/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_6_double_q2_2("ELL  Product double opencl_gpu L6, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q2_2("ELL  Product double opencl_gpu L7, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q2_2("ELL  Product double opencl_gpu L8, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q2_2("ELL  Product double opencl_gpu L9, q2 sort 2", "testdata/poisson_advanced/q2_sort_2/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_6_double_q2_3("ELL  Product double opencl_gpu L6, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q2_3("ELL  Product double opencl_gpu L7, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q2_3("ELL  Product double opencl_gpu L8, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q2_3("ELL  Product double opencl_gpu L9, q2 sort 3", "testdata/poisson_advanced/q2_sort_3/A_9.ell", 10);

ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_6_double_q2_4("ELL  Product double opencl_gpu L6, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_6.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_7_double_q2_4("ELL  Product double opencl_gpu L7, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_7.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_8_double_q2_4("ELL  Product double opencl_gpu L8, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_8.ell", 10);
ProductELLFileBenchmark<tags::OpenCL::GPU, double> opencl_gpu_pareng_9_double_q2_4("ELL  Product double opencl_gpu L9, q2 sort 4", "testdata/poisson_advanced/q2_sort_4/A_9.ell", 10);
#endif
#endif
