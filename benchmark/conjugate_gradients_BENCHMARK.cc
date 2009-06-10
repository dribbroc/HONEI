/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/matrix_io.hh>
#include <honei/math/conjugate_gradients.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class SMELLCGBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        SMELLCGBench(const std::string & id, unsigned long size, unsigned long count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/5pt_10x10.mtx";
            //std::string filename = "/home/user/mgeveler/nobackup/feat2/Featflow2/area51/renumbenchmark/l10/test_2.mtx";
            unsigned long non_zeros(MatrixIO<io_formats::MTX>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DataType_> data(non_zeros);

            MatrixIO<io_formats::MTX>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::MTX>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrixELL<DataType_> smatrix2(rows, columns, r, c, data);

            DenseVector<DataType_> x(rows, DataType_(1.2345));

            DenseVector<DataType_> rhs(rows, DataType_(0));

            Product<Tag_>::value(rhs, smatrix2, x);

            DenseVector<DataType_> diag_inverted(rows, DataType_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DataType_(1)/data[i];
            }

            DenseVector<DataType_> result(rhs.size(), DataType_(1));

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        (ConjugateGradients<Tag_, JAC>::value(smatrix2, rhs, result, diag_inverted, 10ul));
                        );
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
            }
            BenchmarkInfo info;
            BenchmarkInfo info_pre;
            info_pre.flops = ((2 * non_zeros + 2 * rows));
            info.flops = non_zeros * 2 + 12 * rows + 2;
            evaluate(info * 30 + info_pre);
        }
};
SMELLCGBench<tags::CPU, float> SMELLDVPBench_float("SM ELL CG Benchmark CPU float: " , 0, 10);
#ifdef HONEI_CUDA
SMELLCGBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda("SM ELL CG Benchmark CUDA float: " , 0, 6);
#ifdef HONEI_CUDA_DOUBLE
SMELLCGBench<tags::GPU::CUDA, double> SMELLDVPBench_double_cuda("SM ELL CG Benchmark CUDA double: " , 0, 10);
#endif
#endif
