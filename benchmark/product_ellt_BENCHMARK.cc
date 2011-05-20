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
        unsigned long _threads;
        unsigned long _blocks;

    public:
        ProductELLFileBenchmark(const std::string & tag, std::string filename, unsigned long count, unsigned long threads, unsigned long blocks) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _file_name = filename;
            _count = count;
            _threads = threads;
            _blocks = blocks;
        }

        virtual void run()
        {
            std::string filebase(HONEI_SOURCEDIR);
            filebase += "/honei/math/";
            _file_name = filebase + _file_name;
            SparseMatrixELL<DT_> bsmatrix = MatrixIO<io_formats::ELL>::read_matrix(_file_name, DT_(1));
            SparseMatrix<DT_> bla(bsmatrix);
            SparseMatrixELL<DT_>smatrix(bla, _threads);
            std::cout<<smatrix.num_cols_per_row()<<" "<<smatrix.rows()<<std::endl;
            DenseVector<DT_> x(smatrix.rows());
            DenseVector<DT_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }

            Configuration::instance()->set_value("cuda::product_smell_dv_double", _blocks);

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
            evaluate(info * 100);
            }
            std::cout<<"Non Zero Elements: "<<smatrix.used_elements()<<std::endl;
        }
};
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_0("ELL  Product double cuda L5, q2 sort 2, 1 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 1, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_1("ELL  Product double cuda L5, q2 sort 2, 1 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 1, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_2("ELL  Product double cuda L5, q2 sort 2, 1 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 1, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_21("ELL  Product double cuda L5, q2 sort 2, 1 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 1, 1024);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_3("ELL  Product double cuda L5, q2 sort 2, 2 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 2, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_4("ELL  Product double cuda L5, q2 sort 2, 2 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 2, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_5("ELL  Product double cuda L5, q2 sort 2, 2 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 2, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_51("ELL  Product double cuda L5, q2 sort 2, 2 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 2, 1024);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_6("ELL  Product double cuda L5, q2 sort 2, 4 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 4, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_7("ELL  Product double cuda L5, q2 sort 2, 4 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 4, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_8("ELL  Product double cuda L5, q2 sort 2, 4 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 4, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_81("ELL  Product double cuda L5, q2 sort 2, 4 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 4, 1024);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_9("ELL  Product double cuda L5, q2 sort 2, 8 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 8, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_10("ELL  Product double cuda L5, q2 sort 2, 8 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 8, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_11("ELL  Product double cuda L5, q2 sort 2, 8 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 8, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_111("ELL  Product double cuda L5, q2 sort 2, 8 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 8, 1024);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_12("ELL  Product double cuda L5, q2 sort 2, 16 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 16, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_13("ELL  Product double cuda L5, q2 sort 2, 16 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 16, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_14("ELL  Product double cuda L5, q2 sort 2, 16 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 16, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_141("ELL  Product double cuda L5, q2 sort 2, 16 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 16, 1024);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_15("ELL  Product double cuda L5, q2 sort 2, 32 threads, 128 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 32, 128);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_16("ELL  Product double cuda L5, q2 sort 2, 32 threads, 256 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 32, 256);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_17("ELL  Product double cuda L5, q2 sort 2, 32 threads, 512 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 32, 512);
ProductELLFileBenchmark<tags::GPU::CUDA, double> acuda_pareng_7_double_q2_171("ELL  Product double cuda L5, q2 sort 2, 32 threads, 1024 blocks", "testdata/poisson_advanced/q2_sort_2/A_5.ell", 10, 32, 512);
#endif
#endif
