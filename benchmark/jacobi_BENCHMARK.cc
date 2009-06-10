/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/matrix_io.hh>
#include <honei/math/jacobi.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class SMELLJacobiBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        SMELLJacobiBench(const std::string & id, unsigned long size, unsigned long count) :
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
            //filename += "/honei/math/testdata/test_0.mtx";
            unsigned long non_zeros(MatrixIO::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DataType_> data(non_zeros);

            MatrixIO::read_matrix(filename, r, c, data);
            MatrixIO::get_sizes(filename, columns, rows, ax, bx);
            SparseMatrixELL<DataType_> smatrix2(rows, columns, r, c, data);

            DenseVector<DataType_> x(rows, DataType_(1.2345));

            DenseVector<DataType_> rhs(rows, DataType_(0));

            Product<Tag_>::value(rhs, smatrix2, x);

            DenseVector<DataType_> initial_guess(rows, DataType_(1));
            DenseVector<DataType_> diag_inverted(rows, DataType_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DataType_(1)/data[i];
            }

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        Jacobi<Tag_>::value(initial_guess, smatrix2, rhs, 100ul, DataType_(0.7), diag_inverted);
                        );
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
            }
            BenchmarkInfo info;
            info.flops = non_zeros * 2 + 3 * rows;
            evaluate(info * 100);
        }
};
SMELLJacobiBench<tags::CPU, float> SMELLDVPBench_float("SM ELL Jacobi Benchmark CPU float: " , 0, 10);
#ifdef HONEI_CUDA
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda("SM ELL Jacobi Benchmark CUDA float: " , 0, 10);
#ifdef HONEI_CUDA_DOUBLE
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double_cuda("SM ELL Jacobi Benchmark CUDA double: " , 0, 10);
#endif
#endif
