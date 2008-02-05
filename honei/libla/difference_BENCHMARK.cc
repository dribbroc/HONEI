/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/difference.hh>

using namespace honei;

template <typename Tag_, typename DataType_>
class DenseVectorDifferenceBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseVectorDifferenceBench(const std::string & id, unsigned long size, int count) :
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
            for(int i(0) ; i < _count ; ++i)
            {
                BENCHMARK(Difference<Tag_>::value(dv0, dv1));
            }
            BenchmarkInfo info(Difference<>::get_benchmark_info(dv0, dv1));
            evaluate(info);
        }
};
DenseVectorDifferenceBench<tags::CPU, float> DVDBenchfloat1("Dense Vector Difference Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorDifferenceBench<tags::CPU, double> DVDBenchdouble1("Dense Vector Difference Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
DenseVectorDifferenceBench<tags::CPU::MultiCore, float> MCDVDBenchfloat1("MC Dense Vector Difference Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorDifferenceBench<tags::CPU::MultiCore, double> MCDVDBenchdouble1("MC Dense Vector Difference Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#ifdef HONEI_SSE
DenseVectorDifferenceBench<tags::CPU::SSE, float> SSEDVDBenchfloat1("SSE Dense Vector Difference Benchmark - vector size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorDifferenceBench<tags::CPU::SSE, double> SSEDVDBenchdouble1("SSE Dense Vector Difference Benchmark - vector size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
DenseVectorDifferenceBench<tags::CPU::MultiCore::SSE, float> MCSSEDVDBenchfloat1("MC::SSE Dense Vector Difference Benchmark - vector size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorDifferenceBench<tags::CPU::MultiCore::SSE, double> MCSSEDVDBenchdouble1("MC::SSE Dense Vector Difference Benchmark - vector size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
DenseVectorDifferenceBench<tags::Cell, float> CELLDVDBenchfloat1("Cell Dense Vector Difference Benchmark - vector size: 64^4, float", 64ul*64*64*64, 10);
DenseVectorDifferenceBench<tags::Cell, double> CELLDVDBenchdouble1("Cell Dense Vector Difference Benchmark - vector size: 64^4, double", 64ul*64*64*64, 10);
#endif

template <typename Tag_, typename DataType_>
class BandedMatrixDenseMatrixDifferenceBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        BandedMatrixDenseMatrixDifferenceBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size, DataType_(rand()));
            BandedMatrix<DataType_> bm(_size, dv);
            DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
            for(int i(0) ; i < _size ; i += 10)
            {
                bm.insert_band(i, dv);
                bm.insert_band(-i, dv);
            }
            for(int i(0) ; i < _count ; ++i)
            {
                BENCHMARK(Difference<Tag_>::value(bm, dm));
            }
            BenchmarkInfo info(Difference<>::get_benchmark_info(bm, dm));
            evaluate(info);
        }
};
BandedMatrixDenseMatrixDifferenceBench<tags::CPU, float> BMDMDBenchfloat1("Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, float", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU, double> BMDMDBenchdouble1("Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, double", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::MultiCore, float> MCBMDMDBenchfloat1("MC Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, float", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::MultiCore, double> MCBMDMDBenchdouble1("MC Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::SSE, float> SSEBMDMDBenchfloat1("SSE Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, float", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::SSE, double> SSEBMDMDBenchdouble1("SSE Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, double", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::MultiCore::SSE, float> MCSSEBMDMDBenchfloat1("MC::SSE Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, float", 4096, 10);
BandedMatrixDenseMatrixDifferenceBench<tags::CPU::MultiCore::SSE, double> MCSSEBMDMDBenchdouble1("MC::SSE Matrix Difference Benchmark banded/dense - Matrix size: 4096x4096, double", 4096, 10);
#endif
