/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/product.hh>
#include <iostream>

//using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        BandedMatrixDenseVectorProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv1(_size, DataType_(2));
            BandedMatrix<DataType_> bm1(_size, dv1);
            DenseVector<DataType_> dv4(_size, DataType_(3));
            DenseVector<DataType_> dv5(dv4.copy());
            for (int i = 1; i < 14 && i < _size; i++)
            {
                bm1.insert_band(i * 3, dv4.copy());
                bm1.insert_band(-1 * 5 * i, dv5.copy());
            }
            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(Product<Tag_>::value(bm1, dv2));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info);
        }
};

BandedMatrixDenseVectorProductBench<tags::CPU, float> BMDVPBenchfloat("Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64, 10);
BandedMatrixDenseVectorProductBench<tags::CPU, double> BMDVPBenchdouble("Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64*64, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore, float> MCBMDVPBenchfloat("MC Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore, double> MCBMDVPBenchdouble("MC Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64*64, 10);
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, float> SSEBMDVPBenchfloat("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, double> SSEBMDVPBenchdouble("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> MCSSEBMDVPBenchfloat("MC::SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> MCSSEBMDVPBenchdouble("MC::SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
//BandedMatrixDenseVectorProductBench<tags::Cell, float> CELLBMDVPBenchfloat("CELL Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
//BandedMatrixDenseVectorProductBench<tags::Cell, double> CELLBMDVPBenchdouble("CELL Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 10);
#endif


template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductBenchRelax :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        BandedMatrixDenseVectorProductBenchRelax(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv1(_size, DataType_(2));
            BandedMatrix<DataType_> bm1(_size, dv1);
            DenseVector<DataType_> dv4(_size, DataType_(3));
            DenseVector<DataType_> dv5(dv4.copy());
            bm1.insert_band(3, dv4.copy());
            bm1.insert_band(-3, dv5.copy());
            bm1.insert_band(15, dv5.copy());
            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(Product<Tag_>::value(bm1, dv2));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info);
        }
};

#ifdef HONEI_SSE
BandedMatrixDenseVectorProductBenchRelax<tags::CPU::SSE, float> SSEBMDVPBenchfloatRelax("SSE Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBenchRelax<tags::CPU::SSE, double> SSEBMDVPBenchdoubleRelax("SSE Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductBenchRelax<tags::Cell, float> CELLBMDVPBenchfloatRelax("CELL Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBenchRelax<tags::Cell, double> CELLBMDVPBenchdoubleRelax("CELL Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 10);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixProductBench(const std::string & id, unsigned long size, int count) :
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
                BENCHMARK(Product<Tag_>::value(dm0, dm1));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(dm0, dm1));
            evaluate(info);
        }
};
DenseMatrixProductBench<tags::CPU, float> DMPBenchfloat2("Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::CPU, double> DMPBenchdouble2("Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);
DenseMatrixProductBench<tags::CPU::MultiCore, float> MCDMPBenchfloat2("MC Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::CPU::MultiCore, double> MCDMPBenchdouble2("MC Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);
#ifdef HONEI_SSE
DenseMatrixProductBench<tags::CPU::SSE, float> SSEDMPBenchfloat2("SSE Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::CPU::SSE, double> SSEDMPBenchdouble2("SSE Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, float> SSEDMPBenchfloat2SSE("SSE MC Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, double> SSEDMPBenchdouble2SSE("SSE MC Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);
#endif
#ifdef HONEI_CELL
DenseMatrixProductBench<tags::Cell, float> CELLDMPBenchfloat2("CELL Matrix Product Benchmark dense/dense - matrix size: 128x128, float", 128, 10);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixDenseVectorProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
            DenseVector<DataType_> dv0(_size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Product<Tag_>::value(dm0, dv0));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(dm0, dv0));
            evaluate(info);
        }
};
DenseMatrixDenseVectorProductBench<tags::CPU, float> DMDVPBenchfloat("Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, float", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU, double> DMDVPBenchdouble("Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, double", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU::MultiCore, float> DMDVPBenchfloatMC("MC: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, float", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU::MultiCore, double> DMDVPBenchdoubleMC("MC: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, double", 64ul*64, 10);
#ifdef HONEI_SSE
DenseMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> DMDVPBenchfloatMCSSE("MC::SSE: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, float", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> DMDVPBenchdoubleMCSSE("MC::SSE : Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, double", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU::SSE, float> DMDVPBenchfloatSSE("SSE: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, float", 64ul*64, 10);
DenseMatrixDenseVectorProductBench<tags::CPU::SSE, double> DMDVPBenchdoubleSSE("SSE: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, double", 64ul*64, 10);
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SparseMatrixProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            SparseMatrix<DataType_> sm(_size, _size, (unsigned long)((_size*_size)/10));
            for (typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(rand());
                }
            }
            DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(Product<Tag_>::value(sm, dm));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(sm, dm));
            evaluate(info);
        }
};
SparseMatrixProductBench<tags::CPU, float> SMPBenchfloat2("Matrix Product Benchmark sparse/dense - matrix size: 256x256, float", 256, 10);
SparseMatrixProductBench<tags::CPU, double> SMPBenchdouble2("Matrix Product Benchmark sparse/dense - matrix size: 256x256, double", 256, 10);

template <typename Tag_, typename DataType_>
class BandedMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        BandedMatrixProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            DenseVector<DataType_> dv(_size, DataType_(rand()));
            BandedMatrix<DataType_> bm(_size, dv);
            bm.insert_band(1, dv);
            bm.insert_band(-1, dv);
            bm.insert_band(2, dv);
            bm.insert_band(-2, dv);
            bm.insert_band(5, dv);
            bm.insert_band(-5, dv);
            DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            { 
                BENCHMARK(Product<Tag_>::value(bm, dm));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm, dm));
            evaluate(info);
        }
};
BandedMatrixProductBench<tags::CPU, float> BMPBenchfloat2("Matrix Product Benchmark banded/dense - matrix size: 256x256, float", 256, 10);
BandedMatrixProductBench<tags::CPU, double> BMPBenchdouble2("Matrix Product Benchmark banded/dense - matrix size: 256x256, double", 256, 10);
