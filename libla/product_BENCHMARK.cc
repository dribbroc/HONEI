/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/product.hh>
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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv1(_size, DataType_(2));
            BandedMatrix<DataType_> bm1(_size, dv1);
            DenseVector<DataType_> dv4(_size, DataType_(3));
            DenseVector<DataType_> dv5(dv4.copy());
            bm1.insert_band(1, dv4);
            bm1.insert_band(-1, dv5);
            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(Product<Tag_>::value(bm1, dv2));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info<BandedMatrix<DataType_>, DenseVector<DataType_> >(_size, _size));
            evaluate(info);
        }
};

BandedMatrixDenseVectorProductBench<tags::CPU, float> BMDVPBenchfloat("Banded Matrix Dense Vector Product Benchmark - matrix size: 64^4, float", 64ul*64ul*64ul, 1);
BandedMatrixDenseVectorProductBench<tags::CPU, double> BMDVPBenchdouble("Banded Matrix Dense Vector Product Benchmark - matrix size: 32x32, double", 32, 10);
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, float> SSEBMDVPBenchfloat("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, double> SSEBMDVPBenchdouble("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#endif

template <typename DataType_>
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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            {
                DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
                DenseMatrix<DataType_> dm1(_size, _size, DataType_(rand()));
                BENCHMARK(Product<>::value(dm0, dm1));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info<DenseMatrix<DataType_>, DenseMatrix<DataType_> >(_size, _size, _size));
            evaluate(info);
        }
};
DenseMatrixProductBench<float> DMPBenchfloat2("Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<double> DMPBenchdouble2("Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);

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
            BenchmarkInfo info(Product<>::get_benchmark_info<DenseMatrix<DataType_>, DenseVector<DataType_> >(_size, _size));
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

template <typename DataType_>
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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            { 
                SparseMatrix<DataType_> sm(_size, _size, (unsigned long)((_size*_size)/10));
                for (typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(rand());
                    }
                }                
                DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
                BENCHMARK(Product<>::value(sm, dm));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info<SparseMatrix<DataType_>, DenseMatrix<DataType_> >(_size, _size, _size, (double)0.1));
            evaluate(info);
        }
};
SparseMatrixProductBench<float> SMPBenchfloat2("Matrix Product Benchmark sparse/dense - matrix size: 256x256, float", 256, 10);
SparseMatrixProductBench<double> SMPBenchdouble2("Matrix Product Benchmark sparse/dense - matrix size: 256x256, double", 256, 10);


template <typename DataType_>
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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            { 
                DenseVector<DataType_> dv(_size, DataType_(rand()));
                BandedMatrix<DataType_> bm(_size, dv);
                bm.insert_band(1, dv);
                bm.insert_band(-1, dv);
                bm.insert_band(2, dv);
                bm.insert_band(-2, dv);
                bm.insert_band(5, dv);
                bm.insert_band(-5, dv);
                DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
                BENCHMARK(Product<>::value(bm, dm));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info<BandedMatrix<DataType_>, DenseMatrix<DataType_> >(_size, _size, _size, (double)_size * 7 / (_size * _size)));
            evaluate(info);
        }
};
BandedMatrixProductBench<float> BMPBenchfloat2("Matrix Product Benchmark banded/dense - matrix size: 256x256, float", 256, 10);
BandedMatrixProductBench<double> BMPBenchdouble2("Matrix Product Benchmark banded/dense - matrix size: 256x256, double", 256, 10);

