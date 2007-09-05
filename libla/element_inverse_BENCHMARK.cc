/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/element_inverse.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class DenseMatrixElementInverseBench :
    public Benchmark
{
    private:
        int _size;
        int _count;

    public:
        DenseMatrixElementInverseBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
        {
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(12));
            BENCHMARK(ElementInverse<>::value(dm0));
        }
        evaluate(_size * _size);
    }
};

DenseMatrixElementInverseBench<float>  MEIBenchfloat ("Matrix Element Inverse Benchmark: size: 100x100, float",  100, 10);
DenseMatrixElementInverseBench<double> MEIBenchdouble("Matrix Element Inverse Benchmark: size: 100x100, double", 100, 10);
