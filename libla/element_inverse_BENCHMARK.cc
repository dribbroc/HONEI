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
        unsigned long _size;
        int _count;

    public:
        DenseMatrixElementInverseBench(const std::string & id, unsigned long size, int count) :
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
        BenchmarkInfo info(ElementInverse<>::get_benchmark_info<DataType_>(_size, _size));
        evaluate(info);
    }
};

DenseMatrixElementInverseBench<float>  MEIBenchfloat ("Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double> MEIBenchdouble("Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
