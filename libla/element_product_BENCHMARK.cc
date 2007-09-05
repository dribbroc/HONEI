/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/element_product.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class DenseMatrixElementProductBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseMatrixElementProductBench(const std::string & id, int size, int count) :
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
                BENCHMARK(ElementProduct<DataType_>::value(dm0, dm1));
            }
            evaluate(_size*_size);
        }
};

DenseMatrixElementProductBench<float> DMEPBenchfloat1("Dense Matrix Elementwise Product Benchmark - matrix size: 128x128, float", 128, 10);
DenseMatrixElementProductBench<double> DMEPBenchdouble1("Dense Matrix Elementwise Product Benchmark - matrix size: 128x128, double", 128, 10);
//DenseMatrixElementProductBench<float> DMEPBenchfloat2("Dense Matrix Elementwise Product Benchmark - matrix size: 4096x4096, float", 4096, 10);
//DenseMatrixElementProductBench<double> DMEPBenchdouble2("Dense Matrix Elementwise Product Benchmark - matrix size: 4096x4096, double", 4096, 10);
