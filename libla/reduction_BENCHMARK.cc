/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/reduction.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class DenseMatrixRowSumVectorBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixRowSumVectorBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size);
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
                BENCHMARK(dv = Reduction<rt_sum>::value(dm));
            }
            evaluate(_size*_size);
        }
};
DenseMatrixRowSumVectorBench<float> DMRSVBenchfloat2("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<double> DMRSVBenchdouble2("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
