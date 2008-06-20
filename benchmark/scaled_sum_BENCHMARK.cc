/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/la/scaled_sum.hh>


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
                MemoryArbiter::instance()->write<tags::CPU>(dv0.memid(), dv0.address(), dv0.size() * sizeof(DataType_));
                MemoryArbiter::instance()->release_write<tags::CPU>(dv0.memid());
                MemoryArbiter::instance()->write<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
                MemoryArbiter::instance()->release_write<tags::CPU>(dv1.memid());
                for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long l(0) ; l < 10 ; ++l)
                        {
                        ScaledSum<Tag_>::value(dv0, dv1, b);
                        }
                        );
            }
            BenchmarkInfo info(ScaledSum<>::get_benchmark_info(dv0, dv1, b));
            evaluate(info * 10);
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
#ifdef HONEI_CUDA
DenseVectorScaledSumBench<tags::GPU::CUDA, float>
    CUDADVSSBenchfloat1("CUDA Dense Vector ScaledSum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
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
