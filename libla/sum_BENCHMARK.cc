/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/sum.hh>


using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>

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
            register_tag(Tag_::name);
            _size  = size;
            _count = count;
        }



        virtual void run()
        {
            DataType_ alpha = DataType_ (2);
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(42));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Sum<Tag_>::value(dm0, DataType_ (alpha)));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info(dm0, alpha));
            evaluate(info);
    }
};
ScalarDenseMatrixSumBench<tags::CPU, float>  SDMSBenchfloat ("MatrixShift Benchmark - size: 4096x4096, float",  4096, 12);
ScalarDenseMatrixSumBench<tags::CPU, double> SDMSBenchdouble("MatrixShift Benchmark - size: 4096x4096, double", 4096, 12);




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
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv0(_size, DataType_(rand()));
            DenseVector<DataType_> dv1(_size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Sum<Tag_>::value(dv0, dv1));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info(dv0, dv1));
            evaluate(info);
        }
};

DenseVectorSumBench<tags::CPU, float> DVSBenchfloat1("Dense Vector Sum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorSumBench<tags::CPU, double> DVSBenchdouble1("Dense Vector Sum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
DenseVectorSumBench<tags::CPU::SSE, float> SSEDVSBenchfloat1("SSE Dense Vector Sum Benchmark - vector size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorSumBench<tags::CPU::SSE, double> SSEDVSBenchdouble1("SSE Dense Vector Sum Benchmark - vector size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
DenseVectorSumBench<tags::Cell, float> dvs_bench_cell_float("Cell Dense Vector Sum Benchmark - vector size: 64^4, float",
        64ul * 64 * 64 * 64, 10);
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseVectorRangeSumBench(const std::string & id, unsigned long size, int count) :
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
            DenseVectorRange<DataType_> dvr0 (dv0, _size - 2, 0);
            DenseVectorRange<DataType_> dvr1 (dv1, _size - 2, 1);
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Sum<Tag_>::value(dvr0, dvr1));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info(dvr0, dvr1));
            evaluate(info);
        }
};

DenseVectorRangeSumBench<tags::CPU, float> DVRSBenchfloat1("Dense Vector Range Sum Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorRangeSumBench<tags::CPU, double> DVRSBenchdouble1("Dense Vector Range Sum Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
DenseVectorRangeSumBench<tags::CPU::SSE, float> SSEDVRSBenchfloat1("SSE Dense Vector Range Sum Benchmark - vector size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorRangeSumBench<tags::CPU::SSE, double> SSEDVRSBenchdouble1("SSE Dense Vector Range Sum Benchmark - vector size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixSumBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
            DenseMatrix<DataType_> dm1(_size, _size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Sum<Tag_>::value(dm0, dm1));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info(dm0, dm1));
            evaluate(info);
        }
};

DenseMatrixSumBench<tags::CPU, float> DMSBenchfloat1("Dense Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixSumBench<tags::CPU, double> DMSBenchdouble1("Dense Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);
DenseMatrixSumBench<tags::CPU::MultiCore, float> DMSBenchfloat1mc("MC: Dense Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixSumBench<tags::CPU::MultiCore, double> DMSBenchdouble1mc("MC: Dense Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);
#ifdef HONEI_SSE
DenseMatrixSumBench<tags::CPU::SSE, float> DMSBenchfloat1sse("SSE Dense Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixSumBench<tags::CPU::SSE, double> DMSBenchdouble1sse("SSE Dense Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);
DenseMatrixSumBench<tags::CPU::MultiCore::SSE, float> DMSBenchfloat1mcsse("MC SSE : Dense Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixSumBench<tags::CPU::MultiCore::SSE, double> DMSBenchdouble1mcsse("MC SSE: Dense Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);
#endif
#ifdef HONEI_CELL
DenseMatrixSumBench<tags::Cell, float> DMSBenchfloat1cell("Cell Dense Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixBandedMatrixSumBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixBandedMatrixSumBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
            DenseVector<DataType_> dv (_size, DataType_(rand()));
            BandedMatrix<DataType_> bm(_size, dv);
            for(int i = 0 ; i < _size; i += 10) {
                bm.insert_band(i, dv);
                bm.insert_band(-i, dv);
            }
            
            for(int i = 0 ; i < _count ; ++i)
            {
                BENCHMARK(Sum<Tag_>::value(dm0, bm));
            }
            BenchmarkInfo info(Sum<>::get_benchmark_info(dm0, bm));
            evaluate(info);
        }
};

DenseMatrixBandedMatrixSumBench<tags::CPU, float> DMBMSBenchfloat1("Dense Matrix Banded Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixBandedMatrixSumBench<tags::CPU, double> DMBMSBenchdouble1("Dense Matrix Banded Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);
DenseMatrixBandedMatrixSumBench<tags::CPU::MultiCore, float> DMBMSBenchfloat1mc("MC: Dense Matrix Banded Matrix Sum Benchmark - Matrix size: 4048x4048, float", 4048, 10);
DenseMatrixBandedMatrixSumBench<tags::CPU::MultiCore, double> DMSBMBenchdouble1mc("MC: Dense Matrix Banded Matrix Sum Benchmark - Matrix size: 4048x4048, double", 4048, 10);



template <typename DT_, typename Tag_>
class DenseVectorSumBenchTestPlot :
    public Benchmark
{
    public:
        DenseVectorSumBenchTestPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<int> cores;
            for (unsigned long j(1) ; j < 11 ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv0(j * 1000000, DT_(rand()));
                DenseVector<DT_> dv1(j * 1000000, DT_(rand()));
                for(int i = 0; i < 10; ++i)
                {
                    BENCHMARK(Sum<Tag_>::value(dv0, dv1));
                }
                info = Sum<>::get_benchmark_info(dv0, dv1);
                infolist.push_back(info);  
            }
            evaluate_to_plotfile(infolist, cores, 10);
        }
};
DenseVectorSumBenchTestPlot<double, tags::CPU> DVSBTP("Dense Vector Sum Benchmark - vector size: x * 100.000, 0<x<11 - double");
DenseVectorSumBenchTestPlot<double, tags::CPU::MultiCore> MCDVSBTP("MC: Dense Vector Sum Benchmark - vector size: x * 100.000, 0<x<11 - double");
