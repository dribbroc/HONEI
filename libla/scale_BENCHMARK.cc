/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scale.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class ScalarDenseMatrixScaleBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;

    public:
        ScalarDenseMatrixScaleBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = 2.0;
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(23));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Scale<>::value(DataType_ (alpha), dm0));
            }
            BenchmarkInfo info(Scale<>::get_benchmark_info(alpha, dm0));
            evaluate(info);
    }
};
ScalarDenseMatrixScaleBench<float>  SMPBenchfloat ("DenseMatrix Scale Benchmark: size: 4096x4096, float",  4096, 10);
ScalarDenseMatrixScaleBench<double> SMPBenchdouble("DenseMatrix Scale Benchmark: size: 4096x4096, double", 4096, 10);
