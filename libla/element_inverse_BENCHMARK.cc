/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/element_inverse.hh>

using namespace std;
using namespace honei;


template <typename DataType_, typename Tag_>

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
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
        }

        virtual void run()
        {
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(ElementInverse<Tag_>::value(dm0));
            }
            BenchmarkInfo info(ElementInverse<>::get_benchmark_info(dm0));
            evaluate(info);
        }
};

DenseMatrixElementInverseBench<float, tags::CPU>  MEIBenchfloat ("Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double, tags::CPU> MEIBenchdouble("Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
DenseMatrixElementInverseBench<float, tags::CPU::MultiCore>  MCMEIBenchfloat ("MC: Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double, tags::CPU::MultiCore> MCMEIBenchdouble("MC: Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
DenseMatrixElementInverseBench<float, tags::CPU::MultiCore::SSE>  MCSSEMEIBenchfloat ("MC SSE Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double, tags::CPU::MultiCore::SSE> MCSSEMEIBenchdouble("MC SSE Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
DenseMatrixElementInverseBench<float, tags::CPU::SSE>  SSEMEIBenchfloat ("SSE Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double, tags::CPU::SSE> SSEMEIBenchdouble("SSE Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
#endif
#ifdef HONEI_CELL
DenseMatrixElementInverseBench<float, tags::Cell>  CELLMEIBenchfloat ("Cell Matrix Element Inverse Benchmark: size: 4096x4096, float",  4096, 10);
DenseMatrixElementInverseBench<double, tags::Cell> CELLMEIBenchdouble("Cell Matrix Element Inverse Benchmark: size: 4096x4096, double", 4096, 10);
#endif
