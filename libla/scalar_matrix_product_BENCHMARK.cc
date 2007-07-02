/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scalar_matrix_product.hh>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class ScalarDenseMatrixProductBench :
    public Benchmark
{
    private:
        int _size;
        int _count;

    public:
        ScalarDenseMatrixProductBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = 2.0;
            for(int i = 0; i < _count; ++i)
            {
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm0(new DenseMatrix<DataType_>(_size, _size, DataType_(23)) );
                BENCHMARK(ScalarMatrixProduct<DataType_>::value(DataType_ (alpha), *dm0));
            }
        evaluate(_size * _size);
    }
};

ScalarDenseMatrixProductBench<float>  SMPBenchfloat ("Matrixscalierung Benchmark: size: 4000x4000, float",  4000, 10);
ScalarDenseMatrixProductBench<double> SMPBenchdouble("Matrixscalierung Benchmark: size: 4000x4000, double", 4000, 10);
