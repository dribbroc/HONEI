/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

#include <string>
#endif

#include <honei/la/scaled_sum.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/gpu_pool.hh>

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
                        for (unsigned long l(0) ; l < 20 ; ++l)
                        {
                        ScaledSum<Tag_>::value(dv0, dv1, b);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(ScaledSum<>::get_benchmark_info(dv0, dv1, b));
            evaluate(info * 20);
        }
};

DenseVectorScaledSumBench<tags::CPU, float> DVSSBenchfloat1("Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU, double> DVSSBenchdouble1("Dense Vector ScaledSum Benchmark - vector size: 10,000, double", 10000, 10);
DenseVectorScaledSumBench<tags::CPU::MultiCore, float>
    MCDVSSBenchfloat1("MC Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU::MultiCore, double>
    MCDVSSBenchdouble1("MC Dense Vector ScaledSum Benchmark - vector size: 10,000, double", 10000, 10);
#ifdef HONEI_SSE
DenseVectorScaledSumBench<tags::CPU::SSE, float>
    SSEDVSSBenchfloat1("SSE Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU::SSE, double>
    SSEDVSSBenchdouble1("SSE Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU::MultiCore::SSE, float>
    MCSSEDVSSBenchfloat1("MC SSE Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU::MultiCore::SSE, double>
    MCSSEDVSSBenchdouble1("MC SSE Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif
#ifdef HONEI_ITANIUM
DenseVectorScaledSumBench<tags::CPU::Itanium, float>
    ITANIUMDVSSBenchfloat1("Itanium Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::CPU::Itanium, double>
    ITANIUMDVSSBenchdouble1("Itanium Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif
#ifdef HONEI_CUDA
DenseVectorScaledSumBench<tags::GPU::CUDA, float>
    CUDADVSSBenchfloat1("CUDA Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::GPU::MultiCore::CUDA, float>
    MCCUDADVSSBenchfloat1("MC CUDA Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumBench<tags::GPU::CUDA, double>
    CUDADVSSBenchdouble("CUDA Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::GPU::MultiCore::CUDA, double>
    MCCUDADVSSBenchdouble("MC CUDA Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorScaledSumBench<tags::OpenCL::CPU, float>
    OCLCPUDVSSBenchfloat("OpenCL CPU Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::OpenCL::CPU, double>
    OCLCPUDVSSBenchdouble("OpenCL CPU Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::OpenCL::GPU, float>
    OCLGPUDVSSBenchfloat("OpenCL GPU Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumBench<tags::OpenCL::GPU, double>
    OCLGPUDVSSBenchdouble("OpenCL GPU Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif
#endif
#ifdef HONEI_CELL
DenseVectorScaledSumBench<tags::Cell, float> CellDVSSBenchfloat1("CELL Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorScaledSumBench<tags::Cell, double> CellDVSSBenchdouble1("CELL Dense Vector ScaledSum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif

template <typename DT_, typename Tag_>
class DenseVectorScaledSumSPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseVectorScaledSumSPUPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            int temp(Configuration::instance()->get_value("cell::scaled_sum_dense_dense_float", 4));
            int temp2(Configuration::instance()->get_value("cell::scaled_sum_dense_dense_double", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < 81 ; k+=5)
                {
                    Configuration::instance()->set_value("cell::scaled_sum_dense_dense_float", j);
                    Configuration::instance()->set_value("cell::scaled_sum_dense_dense_double", j);
                    cores.push_back(stringify(j) +"SPUs" );
                    DenseVector<DT_> dv0(k * 150000, DT_(rand()));
                    DenseVector<DT_> dv1(k * 150000, DT_(rand()));
                    DT_ alpha(rand());

                    for(int i(0) ; i < 20 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                ScaledSum<Tag_>::value(dv0, dv1, alpha);
                                }
                                );
                    }
                    info = ScaledSum<>::get_benchmark_info(dv0, dv1, alpha);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::scaled_sum_dense_dense_float", temp);
            Configuration::instance()->set_value("cell::scaled_sum_dense_dense_double", temp2);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifdef HONEI_CELL
DenseVectorScaledSumSPUPlot<float, tags::Cell> DVSSSPUF("Cell Dense Vector ScaledSum Benchmark - SPU Count: 1 to 6 - float");
DenseVectorScaledSumSPUPlot<double, tags::Cell> DVSSSPUD("Cell Dense Vector ScaledSum Benchmark - SPU Count: 1 to 6 - double");
#endif

template <typename DT_>
class DenseVectorScaledSumVSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseVectorScaledSumVSPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(tags::CPU::SSE::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            // mc::sse
            for (unsigned long j(1) ; j < 95 ; j+=5)
            {
                cores.push_back(tags::CPU::MultiCore::SSE::name);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                DT_ alpha(rand());

                for(int i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK((ScaledSum<tags::CPU::MultiCore::SSE>::value(dv0, dv1, alpha)));
                }
                info = ScaledSum<>::get_benchmark_info(dv0, dv1, alpha);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long j(1) ; j < 95 ; j+=5)
            {
                cores.push_back(tags::CPU::SSE::name);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                DT_ alpha(rand());

                for(int i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK((ScaledSum<tags::CPU::SSE>::value(dv0, dv1, alpha)));
                }
                info = ScaledSum<>::get_benchmark_info(dv0, dv1, alpha);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
DenseVectorScaledSumVSPlot<float> DVSSVSF("MC vs SSE DenseVector ScaledSum Benchmark - float");
DenseVectorScaledSumVSPlot<double> DVSSVSD("MC vs SSE DenseVector ScaledSum Benchmark - double");
#endif
