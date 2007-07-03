/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>
#include <libla/dense_vector.hh>
#include <tr1/memory>
#include <string>
#endif

#include <libla/matrix_row_sum_vector.hh>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class DenseMatrixRowSumVectorBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseMatrixRowSumVectorBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(_size));
            for(int i = 0; i < _count; ++i)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>(_size, _size, DataType_(rand())) );
                BENCHMARK(*dv = MatrixRowSumVector<DataType_>::value(*dm));
            }
            evaluate(_size*_size);
        }
};

DenseMatrixRowSumVectorBench<float> DMRSVBenchfloat1("Dense Matrix Row Sum Vector Benchmark - matrix size: 2048x2048, float", 2048, 10);
DenseMatrixRowSumVectorBench<double> DMRSVBenchdouble1("Dense Matrix Row Sum Vector Benchmark - matrix size: 2048x2048, double", 2048, 10);
DenseMatrixRowSumVectorBench<float> DMRSVBenchfloat2("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<double> DMRSVBenchdouble2("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
