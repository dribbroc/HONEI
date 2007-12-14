/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <string>
#include <tr1/memory>
#endif

#include <libla/dot_product.hh>

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
                BENCHMARK(p0 = DotProduct<Tag_>::value(dv0, dv1));
            }
            BenchmarkInfo info(DotProduct<>::get_benchmark_info(dv1, dv0));
            evaluate(info);
        }
};
DotProductBench<float, tags::CPU> DPBenchfloat("Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU> DPBenchdouble("Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
DotProductBench<float, tags::CPU::MultiCore> MCDPBenchfloat("MC: Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::MultiCore> MCDPBenchdouble("MC: Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
DotProductBench<float, tags::CPU::MultiCore::SSE> SSEMCDPBenchfloat("MC: SSE Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::MultiCore::SSE> SSEMCDPBenchdouble("MC: SSE  Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
DotProductBench<float, tags::CPU::SSE> SSEDPBenchfloat("SSE Dot Product Benchmark dense/dense - vector size: 64^4 float", 64ul*64*64*64, 10);
DotProductBench<double, tags::CPU::SSE> SSEDPBenchdouble("SSE Dot Product Benchmark dense/dense - vector size: 64^4 double", 64ul*64*64*64, 10);
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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            SparseVector<DataType_> sv(_size, (unsigned long)(_size/10));
            for (typename Vector<DataType_>::ElementIterator i(sv.begin_elements()), i_end(sv.end_elements()) ; i != i_end ; ++i) 
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
SparseDotProductBench<float, tags::CPU::MultiCore> MCSDPBenchfloat("MC: Dot Product Benchmark sparse/dense - vector size: 10,000,000 float", 10000000, 10);
SparseDotProductBench<double, tags::CPU::MultiCore> MCSDPBenchdouble("MC: Dot Product Benchmark sparse/dense - vector size: 10,000,000 double", 10000000, 10);

