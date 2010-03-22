/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <string>
#include <tr1/memory>
#endif

#include <honei/la/dot_product.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>

using namespace std;
using namespace honei;


template <typename DataType_, typename Tag_>

class DotProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DotProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv0(_size, DataType_(rand()));
            DenseVector<DataType_> dv1(_size, DataType_(rand()));
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                        p0 = DotProduct<Tag_>::value(dv0, dv1);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(DotProduct<>::get_benchmark_info(dv1, dv0));
            evaluate(info * 10);
        }
};
DotProductBench<float, tags::CPU> DPBenchfloat("Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU> DPBenchdouble("Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
DotProductBench<float, tags::CPU::MultiCore> MCDPBenchfloat("MC: Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::MultiCore> MCDPBenchdouble("MC: Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
DotProductBench<float, tags::CPU::SSE> SSEDPBenchfloat("SSE Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::SSE> SSEDPBenchdouble("SSE Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
DotProductBench<float, tags::CPU::MultiCore::SSE> SSEMCDPBenchfloat("MC: SSE Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::MultiCore::SSE> SSEMCDPBenchdouble("MC: SSE  Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
#endif
#ifdef HONEI_CUDA
DotProductBench<float, tags::GPU::CUDA> CUDADPBenchfloat("CUDA Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
#ifdef HONEI_CUDA_DOUBLE
DotProductBench<double, tags::GPU::CUDA> CUDADPBenchdouble("CUDA Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
#endif
#endif
#ifdef HONEI_CELL
DotProductBench<float, tags::Cell> CELLDPBenchfloat("Cell Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
#endif

template <typename DataType_, typename Tag_>

class SparseDotProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SparseDotProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            SparseVector<DataType_> sv(_size, (unsigned long)(_size/10));
            for (typename SparseVector<DataType_>::ElementIterator i(sv.begin_elements()), i_end(sv.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(rand());
                }
            }
            DenseVector<DataType_> dv(_size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(p0 = DotProduct<Tag_>::value(sv,dv));
            }
            BenchmarkInfo info(DotProduct<>::get_benchmark_info(sv, dv));
            evaluate(info);
        }
};
SparseDotProductBench<float, tags::CPU> SDPBenchfloat("Dot Product Benchmark sparse/dense - vector size: 10,000,000 float", 10000000, 10);
SparseDotProductBench<double, tags::CPU> SDPBenchdouble("Dot Product Benchmark sparse/dense - vector size: 10,000,000 double", 10000000, 10);
//SparseDotProductBench<float, tags::CPU::MultiCore> MCSDPBenchfloat("MC: Dot Product Benchmark sparse/dense - vector size: 10,000,000 float", 10000000, 10);
//SparseDotProductBench<double, tags::CPU::MultiCore> MCSDPBenchdouble("MC: Dot Product Benchmark sparse/dense - vector size: 10,000,000 double", 10000000, 10);


template <typename DT_, typename Tag_>
class DenseVectorDotProductSPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseVectorDotProductSPUPlot(const std::string & id) :
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

            int temp(Configuration::instance()->get_value("cell::dot_product_dense_dense_float", 4));
            int temp2(Configuration::instance()->get_value("cell::dot_product_dense_dense_double", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < 81 ; k+=5)
                {
                    Configuration::instance()->set_value("cell::dot_product_dense_dense_float", j);
                    Configuration::instance()->set_value("cell::dot_product_dense_dense_double", j);
                    cores.push_back(stringify(j) +"SPUs" );
                    DenseVector<DT_> dv0(k * 150000, DT_(rand()));
                    DenseVector<DT_> dv1(k * 150000, DT_(rand()));

                    for(int i(0) ; i < 20 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                DotProduct<Tag_>::value(dv0, dv1);
                                }
                                );
                    }
                    info = DotProduct<>::get_benchmark_info(dv0, dv1);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::dot_product_dense_dense_float", temp);
            Configuration::instance()->set_value("cell::dot_product_dense_dense_double", temp2);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifdef HONEI_CELL
DenseVectorDotProductSPUPlot<float, tags::Cell> DVDPSPUF("Cell Dense Vector DotProduct Benchmark - SPU Count: 1 to 6 - float");
DenseVectorDotProductSPUPlot<double, tags::Cell> DVDPSPUD("Cell Dense Vector DotProduct Benchmark - SPU Count: 1 to 6 - double");
#endif
