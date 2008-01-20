/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/element_product.hh>

using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>

class DenseMatrixElementProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixElementProductBench(const std::string & id, unsigned long size, int count) :
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
            DenseMatrix<DataType_> dm1(_size, _size, DataType_(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(ElementProduct<Tag_>::value(dm0, dm1));
            }
            BenchmarkInfo info(ElementProduct<>::get_benchmark_info(dm0, dm1));
            evaluate(info);
        }
};
DenseMatrixElementProductBench<tags::CPU, float> DMEPBenchfloat1("Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementProductBench<tags::CPU, double> DMEPBenchdouble1("Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixElementProductBench<tags::CPU::MultiCore, float> DMEPBenchfloat1mc("MC: Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementProductBench<tags::CPU::MultiCore, double> DMEPBenchdouble1mc("MC: Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
DenseMatrixElementProductBench<tags::CPU::SSE, float> SSEDMEPBenchfloat1("SSE Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementProductBench<tags::CPU::SSE, double> SSEDMEPBenchdouble1("SSE Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixElementProductBench<tags::CPU::MultiCore::SSE, float> SSEMCDMEPBenchfloat1mc("MC: SSE Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementProductBench<tags::CPU::MultiCore::SSE, double> SSEMCDMEPBenchdouble1mc("MC: SSE Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, double", 4096, 10);
#endif


template <typename Tag_, typename DataType_>

class SparseMatrixElementProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SparseMatrixElementProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            SparseMatrix<DataType_> sm(_size, _size, (unsigned long)(_size/10)); 
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
                BENCHMARK(ElementProduct<Tag_>::value(sm, dm));
            }
            BenchmarkInfo info(ElementProduct<>::get_benchmark_info(sm, dm));
            evaluate(info);
        }
};
SparseMatrixElementProductBench<tags::CPU, float> SMEPBenchfloat1("Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, float", 2048, 10);
SparseMatrixElementProductBench<tags::CPU, double> SMEPBenchdouble1("Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, double", 2048, 10);
SparseMatrixElementProductBench<tags::CPU::MultiCore, float> SMEPBenchfloat1MC("MC: Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, float", 2048, 10);
SparseMatrixElementProductBench<tags::CPU::MultiCore, double> SMEPBenchdouble1MC("MC: Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, double", 2048, 10);

template <typename Tag_, typename DataType_>

class BandedMatrixElementProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        BandedMatrixElementProductBench(const std::string & id, unsigned long size, int count) :
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
                BENCHMARK(ElementProduct<Tag_>::value(bm, dm));
            }
        BenchmarkInfo info(ElementProduct<>::get_benchmark_info(bm, dm));
        evaluate(info);
        }
};
BandedMatrixElementProductBench<tags::CPU, float> BMEPBenchfloat1("Matrix Elementwise Product Benchmark banded/dense - matrix size: 2048x2048, float", 2048, 10);
BandedMatrixElementProductBench<tags::CPU, double> BMEPBenchdouble1("Matrix Elementwise Product Benchmark banded/dense - matrix size: 2048x2048, double", 2048, 10);

template <typename Tag_, typename DataType_>

class DenseVectorElementProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseVectorElementProductBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            DenseVector<DataType_> dv1(_size, DataType_(rand()));
            DenseVector<DataType_> dv2(_size, DataType_(rand()));

            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(ElementProduct<Tag_>::value(dv1, dv2));
            }
        BenchmarkInfo info(ElementProduct<>::get_benchmark_info(dv1, dv2));
        evaluate(info);
        }
};
DenseVectorElementProductBench<tags::CPU, float> DVEPBenchfloat1("Vector Elementwise Product Benchmark dense/dense size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU, double> DVEPBenchdouble1("Vector Elementwise Product Benchmark dense/dense size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU::MultiCore, float> mc_DVEPBenchfloat1("MC: Vector Elementwise Product Benchmark dense/dense size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU::MultiCore, double> mc_DVEPBenchdouble1("MC: Vector Elementwise Product Benchmark dense/dense size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#ifdef HONEI_SSE
DenseVectorElementProductBench<tags::CPU::SSE, float> SSEDVEPBenchfloat1("SSE: Vector Elementwise Product Benchmark dense/dense size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU::SSE, double> SSEDVEPBenchdouble1("SSE: Vector Elementwise Product Benchmark dense/dense size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU::MultiCore::SSE, float> MCSSEDVEPBenchfloat1("MC::SSE: Vector Elementwise Product Benchmark dense/dense size: 64^4, float", 64ul*64ul*64ul*64ul, 10);
DenseVectorElementProductBench<tags::CPU::MultiCore::SSE, double> MCSSEDVEPBenchdouble1("MC::SSE: Vector Elementwise Product Benchmark dense/dense size: 64^4, double", 64ul*64ul*64ul*64ul, 10);
#endif
#ifdef HONEI_CELL
DenseVectorElementProductBench<tags::Cell, float>  CELLDVEPBenchfloat ("Cell Vector Element Product Benchmark dense/dense size: 64^4, float",  64ul*64ul*64ul*64ul, 10);
#endif
