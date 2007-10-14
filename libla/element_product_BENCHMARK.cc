/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/element_product.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

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
                BENCHMARK(ElementProduct<DataType_>::value(dm0, dm1));
            }
            evaluate(_size*_size);
        }
};
DenseMatrixElementProductBench<float> DMEPBenchfloat1("Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementProductBench<double> DMEPBenchdouble1("Matrix Elementwise Product Benchmark dense/dense - matrix size: 4096x4096, double", 4096, 10);


template <typename DataType_>

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
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DataType_ p0;
            for(int i = 0; i < _count; ++i)
            {
                SparseMatrix<DataType_> sm(_size, _size, (unsigned long)(_size/10)); 
                for (typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DataType_(rand());
                    }
                }
                DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
                BENCHMARK(ElementProduct<DataType_>::value(sm, dm));
            }
            evaluate((unsigned long)((_size*_size)/10));
        }
};
SparseMatrixElementProductBench<float> SMEPBenchfloat1("Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, float", 2048, 10);
SparseMatrixElementProductBench<double> SMEPBenchdouble1("Matrix Elementwise Product Benchmark sparse/dense - matrix size: 2048x2048, double", 2048, 10);


template <typename DataType_>

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
                BENCHMARK(ElementProduct<DataType_>::value(bm, dm));
            }
            evaluate((unsigned long)(_size*7-16));
        }
};
BandedMatrixElementProductBench<float> BMEPBenchfloat1("Matrix Elementwise Product Benchmark banded/dense - matrix size: 2048x2048, float", 2048, 10);
BandedMatrixElementProductBench<double> BMEPBenchdouble1("Matrix Elementwise Product Benchmark banded/dense - matrix size: 2048x2048, double", 2048, 10);

