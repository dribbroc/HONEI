/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/product.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class DenseMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
                DenseMatrix<DataType_> dm1(_size, _size, DataType_(rand()));
                BENCHMARK(Product<DataType_>::value(dm0, dm1));
            }
            evaluate(_size*_size*_size*2);
        }
};

DenseMatrixProductBench<float> DMPBenchfloat("Dense Matrix Product Benchmark - matrix size: 32x32, float", 32, 10);
DenseMatrixProductBench<double> DMPBenchdouble("Dense Matrix Product Benchmark - matrix size: 32x32, double", 32, 10);
//DenseMatrixProductBench<float> DMPBenchfloat2("Dense Matrix Product Benchmark - matrix size: 256x256, float", 256, 10);
//DenseMatrixProductBench<double> DMPBenchdouble2("Dense Matrix Product Benchmark - matrix size: 256x256, double", 256, 10);
