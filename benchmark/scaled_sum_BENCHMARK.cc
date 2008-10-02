/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/la/scaled_sum.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>


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
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
                DenseVector<DataType_> dv0(_size, static_cast<DataType_>(rand()));
                DenseVector<DataType_> dv1(_size, static_cast<DataType_>(rand()));
                DataType_ b(1234.56789);
                for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long l(0) ; l < 10 ; ++l)
                        {
                        ScaledSum<Tag_>::value(dv0, dv1, b);
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
                        }
                        );
            }
            BenchmarkInfo info(ScaledSum<>::get_benchmark_info(dv0, dv1, b));
            evaluate(info * 10);
        }
};

DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat1("CPU Dense Vector ScaledSum Benchmark - vector size: 65536, float", 65536ul, 10);
DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat1b("CPU Dense Vector ScaledSum Benchmark - vector size: 655360, float", 655360ul, 10);
DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat2("CPU Dense Vector ScaledSum Benchmark - vector size: 1376256, float", 1376256ul, 10);
DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat3("CPU Dense Vector ScaledSum Benchmark - vector size: 268976, float", 2686976ul, 10);
DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat4("CPU Dense Vector ScaledSum Benchmark - vector size: 3997696, float", 3997696ul, 10);

#ifdef HONEI_SSE
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat1("SSE Dense Vector ScaledSum Benchmark - vector size: 65536, float", 65536ul, 10);
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat1b("SSE Dense Vector ScaledSum Benchmark - vector size: 655360, float", 655360ul, 10);
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat2("SSE Dense Vector ScaledSum Benchmark - vector size: 1376256, float", 1376256ul, 10);
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat3("SSE Dense Vector ScaledSum Benchmark - vector size: 268976, float", 2686976ul, 10);
DenseVectorScaledSumBench<tags::CPU::SSE, float> SSEDVSSBenchfloat4("SSE Dense Vector ScaledSum Benchmark - vector size: 3997696, float", 3997696ul, 10);
#endif

#ifdef HONEI_CUDA
DenseVectorScaledSumBench<tags::GPU::CUDA, float> CUDADVSSBenchfloat1("CUDA Dense Vector ScaledSum Benchmark - vector size: 65536, float", 65536ul, 10);
DenseVectorScaledSumBench<tags::GPU::CUDA, float> CUDADVSSBenchfloat1b("CUDA Dense Vector ScaledSum Benchmark - vector size: 655360, float", 655360ul, 10);
DenseVectorScaledSumBench<tags::GPU::CUDA, float> CUDADVSSBenchfloat2("CUDA Dense Vector ScaledSum Benchmark - vector size: 1376256, float", 1376256ul, 10);
DenseVectorScaledSumBench<tags::GPU::CUDA, float> CUDADVSSBenchfloat3("CUDA Dense Vector ScaledSum Benchmark - vector size: 268976, float", 2686976ul, 10);
DenseVectorScaledSumBench<tags::GPU::CUDA, float> CUDADVSSBenchfloat4("CUDA Dense Vector ScaledSum Benchmark - vector size: 3997696, float", 3997696ul, 10);
#endif

#ifdef HONEI_CELL
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat1("Cell Dense Vector ScaledSum Benchmark - vector size: 65536, float", 65536ul, 10);
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat1b("Cell Dense Vector ScaledSum Benchmark - vector size: 655360, float", 655360ul, 10);
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat2("Cell Dense Vector ScaledSum Benchmark - vector size: 1376256, float", 1376256ul, 10);
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat3("Cell Dense Vector ScaledSum Benchmark - vector size: 268976, float", 2686976ul, 10);
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat4("Cell Dense Vector ScaledSum Benchmark - vector size: 3997696, float", 3997696ul, 10);
#endif
