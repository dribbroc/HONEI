/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/reduction.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
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
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size);
            DenseMatrix<DataType_> dm(_size, _size, static_cast<DataType_>(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(dv = (Reduction<rt_sum, Tag_>::value(dm)));
            }
            BenchmarkInfo info(Reduction<rt_sum>::get_benchmark_info(dm));
            evaluate(info);
        }
};
DenseMatrixRowSumVectorBench<tags::CPU, float> DMRSVBenchfloat1("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU, double> DMRSVBenchdouble1("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore, float> DMRSVBenchfloat2("MC: Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore, double> DMRSVBenchdouble2("MC: Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
