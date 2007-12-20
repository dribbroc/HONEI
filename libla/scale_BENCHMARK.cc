/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scale.hh>

using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>

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
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = 2.0;
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(23));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Scale<Tag_>::value(DataType_ (alpha), dm0));
            }
            BenchmarkInfo info(Scale<>::get_benchmark_info(alpha, dm0));
            evaluate(info);
    }
};
ScalarDenseMatrixScaleBench<tags::CPU, float>             SMPBenchfloat   ("DenseMatrix Scale Benchmark: size: 8192x8192, float",      8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU, double>            SMPBenchdouble  ("DenseMatrix Scale Benchmark: size: 8192x8192, double",     8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore, float>  SMPBenchfloatMC ("MC: DenseMatrix Scale Benchmark: size: 8192x8192, float",  8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore, double> SMPBenchdoubleMC("MC: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
#ifdef HONEI_SSE
ScalarDenseMatrixScaleBench<tags::CPU::SSE, float>
        SMPBenchfloatSSE("SSE: DenseMatrix Scale Benchmark: size: 8192x8192, float", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::SSE, double>
        SMPBenchdoubleSSE("SSE: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore::SSE, float>
        SMPBenchfloatMCSSE("MC SSE: DenseMatrix Scale Benchmark: size: 8192x8192, float", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore::SSE, double>
        SMPBenchdoubleMCSSE("MC SSE: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
#endif
