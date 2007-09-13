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
        int _size;
        int _count;

    public:
        ScalarDenseMatrixSumBench(const std::string & id, int size, int count) :
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
        evaluate(_size * _size);
    }
};

ScalarDenseMatrixSumBench<float>  SDMSBenchfloat ("MatrixShift Benchmark: size: 200x200, float",  200, 12);
ScalarDenseMatrixSumBench<double> SDMSBenchdouble("MatrixShift Benchmark: size: 200x200, double", 200, 12);




template <typename Tag_, typename DataType_>

class DenseVectorSumBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseVectorSumBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
            {
                DenseVector<DataType_> dv0(_size, static_cast<DataType_>(rand()));
                DenseVector<DataType_> dv1(_size, static_cast<DataType_>(rand()));
                BENCHMARK(DenseVector<DataType_> sum1(Sum<Tag_>::value(dv0, dv1)));
            }
            evaluate(3*_size);
        }
};

DenseVectorSumBench<tags::CPU, float> DVSBenchfloat1("Dense Vector Sum Benchmark - vector size: 10,000, float", 10000, 10);
DenseVectorSumBench<tags::CPU, double> DVSBenchdouble1("Dense Vector Sum Benchmark - vector size: 10,000, double", 10000, 10);
#ifdef HONEI_SSE
DenseVectorSumBench<tags::CPU::SSE, float> SSEDVSBenchfloat1("SSE Dense Vector Sum Benchmark - vector size: 640,000, float", 640000, 1);
#endif
//DenseVectorSumBench<float> DVSBenchfloat2("Dense Vector Scaled Sum Benchmark - vector size: 10,000,000, float", 10000000, 10);
//DenseVectorScaledSumBench<double> DVSSBenchdouble2("Dense Vector Scaled Sum Benchmark - vector size: 10,000,000, double", 10000000, 10);
