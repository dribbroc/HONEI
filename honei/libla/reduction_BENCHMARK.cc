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
#ifdef HONEI_SSE
DenseMatrixRowSumVectorBench<tags::CPU::SSE, float> SSEDMRSVBenchfloat1("SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::SSE, double> SSEDMRSVBenchdouble1("SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, float> SSEDMRSVBenchfloat2("MC::SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, double> SSEDMRSVBenchdouble2("MC::SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixRowSumVectorBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SparseMatrixRowSumVectorBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size);
            SparseMatrix<DataType_> sm(_size, _size, (unsigned long)(_size/10));
            for (typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(rand());
                }
            }
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(dv = (Reduction<rt_sum, Tag_>::value(sm)));
            }
            BenchmarkInfo info(Reduction<rt_sum>::get_benchmark_info(sm));
            evaluate(info);
        }
};
SparseMatrixRowSumVectorBench<tags::CPU, float> SMRSVBenchfloat1("Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU, double> SMRSVBenchdouble1("Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore, float> SMRSVBenchfloat2("MC: Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore, double> SMRSVBenchdouble2("MC: Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
SparseMatrixRowSumVectorBench<tags::CPU::SSE, float> SSESMRSVBenchfloat1("SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::SSE, double> SSESMRSVBenchdouble1("SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, float> SSESMRSVBenchfloat2("MC::SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, double> SSESMRSVBenchdouble2("MC::SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#endif
