/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scaled_sum.hh>


using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>
class DenseVectorScaledSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseVectorScaledSumBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
            {
                DenseVector<DataType_> dv0(_size, static_cast<DataType_>(rand()));
                DenseVector<DataType_> dv1(_size, static_cast<DataType_>(rand()));
                DataType_ b(1234.56789);
                BENCHMARK(DenseVector<DataType_> sum1(ScaledSum<Tag_>::value(dv0, dv1, b)));
            }
            evaluate(3*_size);
        }
};

DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat1("Dense Vector ScaledSum Benchmark - vector size: (3*4)^6, float", 3*4000, 10);
DenseVectorScaledSumBench<tags::CPU, double> DVSSBenchdouble1("Dense Vector ScaledSum Benchmark - vector size: 10,000, double", 10000, 10);
#ifdef HONEI_SSE
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat1("SSE Dense Vector ScaledSum Benchmark - vector size: (3*4)^6, float", 3*4000, 10);
#endif
