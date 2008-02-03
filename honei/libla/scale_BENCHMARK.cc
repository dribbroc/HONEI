/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/scale.hh>

using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>

class ScalarDenseVectorScaleBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;

    public:
        ScalarDenseVectorScaleBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = 2.0;
            DenseVector<DataType_> dv(_size, DataType_(23));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Scale<Tag_>::value(dv, DataType_ (alpha)));
            }
            BenchmarkInfo info(Scale<>::get_benchmark_info(dv, alpha));
            evaluate(info);
    }
};
ScalarDenseVectorScaleBench<tags::CPU, float>             SVPBenchfloat   ("Densevector Scale Benchmark - vector size: 64^4, float",      64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU, double>            SVPBenchdouble  ("Densevector Scale Benchmark - vector size: 64^4, double",     64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU::MultiCore, float>  SVPBenchfloatMC ("MC: Densevector Scale Benchmark - vector size: 64^4, float",  64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU::MultiCore, double> SVPBenchdoubleMC("MC: Densevector Scale Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
ScalarDenseVectorScaleBench<tags::CPU::SSE, float>              SVPBenchfloatSSE("SSE: Densevector Scale Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU::SSE, double>             SVPBenchdoubleSSE("SSE: Densevector Scale Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU::MultiCore::SSE, float>   SVPBenchfloatMCSSE("MC SSE: Densevector Scale Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
ScalarDenseVectorScaleBench<tags::CPU::MultiCore::SSE, double>  SVPBenchdoubleMCSSE("MC SSE: Densevector Scale Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif

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
                BENCHMARK(Scale<Tag_>::value(dm0, DataType_ (alpha)));
            }
            BenchmarkInfo info(Scale<>::get_benchmark_info(dm0, alpha));
            evaluate(info);
    }
};
ScalarDenseMatrixScaleBench<tags::CPU, float>             SMPBenchfloat   ("DenseMatrix Scale Benchmark: size: 8192x8192, float",      8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU, double>            SMPBenchdouble  ("DenseMatrix Scale Benchmark: size: 8192x8192, double",     8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore, float>  SMPBenchfloatMC ("MC: DenseMatrix Scale Benchmark: size: 8192x8192, float",  8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore, double> SMPBenchdoubleMC("MC: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
#ifdef HONEI_SSE
ScalarDenseMatrixScaleBench<tags::CPU::SSE, float>              SMPBenchfloatSSE("SSE: DenseMatrix Scale Benchmark: size: 8192x8192, float", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::SSE, double>             SMPBenchdoubleSSE("SSE: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore::SSE, float>   SMPBenchfloatMCSSE("MC SSE: DenseMatrix Scale Benchmark: size: 8192x8192, float", 8192, 10);
ScalarDenseMatrixScaleBench<tags::CPU::MultiCore::SSE, double>  SMPBenchdoubleMCSSE("MC SSE: DenseMatrix Scale Benchmark: size: 8192x8192, double", 8192, 10);
#endif
