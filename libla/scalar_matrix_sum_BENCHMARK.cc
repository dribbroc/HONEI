/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scalar_matrix_sum.hh>

 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class ScalarDenseMatrixSumBench :
    public Benchmark
{
    private:
        int _size;
        int _count;

    public:
        ScalarDenseMatrixSumBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = DataType_ (2);
            for(int i = 0; i < _count; ++i)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm0(new DenseMatrix<DataType_>(_size, _size, DataType_(42)) );
                BENCHMARK(ScalarMatrixSum<DataType_>::value(DataType_ (alpha), *dm0));
            }
        evaluate(_size * _size);
    }
};

ScalarDenseMatrixSumBench<float>  SDMSBenchfloat ("MatrixShift Benchmark: size: 2000x2000, float",  2000, 12);
ScalarDenseMatrixSumBench<double> SDMSBenchdouble("MatrixShift Benchmark: size: 2000x2000, double", 2000, 12);
