/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/sum.hh>


using namespace std;
using namespace honei;


template <typename DataType_>

class ScalarDenseMatrixSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;

    public:
        ScalarDenseMatrixSumBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = DataType_ (2);
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> dm0(_size, _size, DataType_(42));
                BENCHMARK(Sum<>::value(dm0, DataType_ (alpha)));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info<DenseMatrix<DataType_>, DataType_>(1, _size, _size));
            evaluate(info);
    }
};
ScalarDenseMatrixSumBench<float>  SDMSBenchfloat ("MatrixShift Benchmark - size: 4096x4096, float",  4096, 12);
ScalarDenseMatrixSumBench<double> SDMSBenchdouble("MatrixShift Benchmark - size: 4096x4096, double", 4096, 12);




template <typename Tag_, typename DataType_>

class DenseVectorSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseVectorSumBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv0(_size, DataType_(rand()));
            DenseVector<DataType_> dv1(_size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(DenseVector<DataType_> sum1(Sum<Tag_>::value(dv0, dv1)));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info<DenseVector<DataType_>, DenseVector<DataType_> >(0, _size));
            evaluate(info);
        }
};

DenseVectorSumBench<tags::CPU, float> DVSBenchfloat1("Dense Vector Sum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorSumBench<tags::CPU, double> DVSBenchdouble1("Dense Vector Sum Benchmark - vector size: 10,000,000, double", 10000000, 10);
#ifdef HONEI_SSE
DenseVectorSumBench<tags::CPU::SSE, float> SSEDVSBenchfloat1("SSE Dense Vector Sum Benchmark - vector size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
DenseVectorSumBench<tags::Cell, float> dvs_bench_cell_float("Cell Dense Vector Sum Benchmark - vector size: 64^4, float",
        64ul * 64 * 64 * 64, 10);
#endif
