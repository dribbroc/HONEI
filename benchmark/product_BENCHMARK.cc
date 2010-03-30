/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/matrix_io.hh>
#include <honei/la/product.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class GenericQ1MatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        GenericQ1MatrixDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count) :
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

            bm1.insert_band(- (unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size) + 1, dv4.copy());
            bm1.insert_band(-1, dv4.copy());
            bm1.insert_band(0, dv4.copy());
            bm1.insert_band(1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size)+ 1, dv4.copy());

            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(bm1, dv2);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};
GenericQ1MatrixDenseVectorProductBench<tags::CPU, float> GenericQ1DVPBenchfloat("Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU, double> GenericQ1DVPBenchdouble("Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU::MultiCore, float> GenericMCQ1DVPBenchfloat("MC Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025, float", 1025*1025ul, 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU::MultiCore, double> GenericMCQ1DVPBenchdouble("MC Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 64^3, double", 1025ul*1025, 10);
#ifdef HONEI_SSE
GenericQ1MatrixDenseVectorProductBench<tags::CPU::SSE, float> GenericSSEQ1DVPBenchfloat("SSE Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025, float",1025ul * 1025 , 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU::SSE, double> GenericSSEQ1DVPBenchdouble("SSE Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025 ,double", 1025ul*1025, 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> GenericMCSSEQ1DVPBenchfloat("MC::SSE Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025, float", 1025ul*1025, 10);
GenericQ1MatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> GenericMCSSEQ1DVPBenchdouble("MC::SSE Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025, double", 1025*1025, 10);
#endif
#ifdef HONEI_CUDA
GenericQ1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> GenericCUDAQ1DVPBenchfloat("CUDA Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025*1025, float",1025ul * 1025 , 10);
#endif
#ifdef HONEI_CELL
GenericQ1MatrixDenseVectorProductBench<tags::Cell, float> GenericCELLQ1DVPBenchfloat("CELL Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025^2, float", 1025ul * 1025, 10);
GenericQ1MatrixDenseVectorProductBench<tags::Cell, double> GenericCELLQ1DVPBenchdouble("CELL Banded Matrix (Q1 generic) Dense Vector Product Benchmark - matrix size: 1025^2, double", 1025ul * 1025, 10);
#endif

template <typename Tag_, typename DataType_>
class Q1MatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        Q1MatrixDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count) :
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

            bm1.insert_band(- (unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size) + 1, dv4.copy());
            bm1.insert_band(-1, dv4.copy());
            bm1.insert_band(0, dv4.copy());
            bm1.insert_band(1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size)+ 1, dv4.copy());

            BandedMatrixQ1<DataType_> qm1(bm1);
            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(qm1, dv2);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat("Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, double> Q1DVPBenchdouble("Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::MultiCore, float> Q1DVPBenchMCfloat("Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, MC float", 1025ul*1025, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::MultiCore, double> Q1DVPBenchMCdouble("Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, MC double", 1025ul*1025, 10);
#ifdef HONEI_SSE
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> SSEQ1DVPBenchfloat("SSE Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025, float",1025ul * 1025 , 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, double> SSEQ1DVPBenchdouble("SSE Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025 ,double", 1025ul*1025, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> MCSSEQ1DVPBenchMCfloat("MC::SSE Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, MC SSE float", 1025ul*1025, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> MCSSEQ1DVPBenchMCdouble("MC::SSE Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, MC SSEdouble", 1025ul*1025, 10);
#endif
#ifdef HONEI_CUDA
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> CUDAQ1DVPBenchfloat("CUDA Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025, float",1025ul * 1025 , 10);
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, double> CUDAQ1DVPBenchdouble("CUDA Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025, double",1025ul * 1025 , 10);
#endif
#endif
#ifdef HONEI_CELL
Q1MatrixDenseVectorProductBench<tags::Cell, float> CELLQ1DVPBenchfloat("CELL Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025^2, float", 1025ul * 1025, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, double> CELLQ1DVPBenchdouble("CELL Banded Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025^2, double", 1025ul * 1025, 10);
#endif

template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        BandedMatrixDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count) :
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
            for (long i = 1; i < 14 && i < (signed)_size; i++)
            {
                bm1.insert_band(i * 3, dv4.copy());
                bm1.insert_band(-1 * 5 * i, dv5.copy());
            }

            DenseVector<DataType_> dv2(_size, DataType_(4));
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(bm1, dv2);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};

BandedMatrixDenseVectorProductBench<tags::CPU, float> BMDVPBenchfloat("Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64, 100);
BandedMatrixDenseVectorProductBench<tags::CPU, double> BMDVPBenchdouble("Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64*64, 100);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore, float> MCBMDVPBenchfloat("MC Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 100);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore, double> MCBMDVPBenchdouble("MC Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64*64, 100);
#ifdef HONEI_SSE
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, float> SSEBMDVPBenchfloat("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64, 100);
BandedMatrixDenseVectorProductBench<tags::CPU::SSE, double> SSEBMDVPBenchdouble("SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 100);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> MCSSEBMDVPBenchfloat("MC::SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 100);
BandedMatrixDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> MCSSEBMDVPBenchdouble("MC::SSE Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 100);
#endif
#ifdef HONEI_CUDA
BandedMatrixDenseVectorProductBench<tags::GPU::CUDA, float> CUDABMDVPBenchfloat("CUDA Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64, 10);
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductBench<tags::Cell, float> CELLBMDVPBenchfloat("CELL Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 100);
BandedMatrixDenseVectorProductBench<tags::Cell, double> CELLBMDVPBenchdouble("CELL Banded Matrix Dense Vector Product Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 100);
#endif


template <typename Tag_, typename DataType_>
class BandedMatrixDenseVectorProductBenchRelax :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        BandedMatrixDenseVectorProductBenchRelax(const std::string & id, unsigned long size, unsigned long count) :
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
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(bm1, dv2);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};

#ifdef HONEI_SSE
BandedMatrixDenseVectorProductBenchRelax<tags::CPU::SSE, float> SSEBMDVPBenchfloatRelax("SSE Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, float", 64ul*64*64, 100);
BandedMatrixDenseVectorProductBenchRelax<tags::CPU::SSE, double> SSEBMDVPBenchdoubleRelax("SSE Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 100);
#endif
#ifdef HONEI_CUDA
BandedMatrixDenseVectorProductBenchRelax<tags::GPU::CUDA, float> CUDABMDVPBenchfloatRelax("CUDA Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, float", 64ul*64*64, 10);
#endif
#ifdef HONEI_CELL
BandedMatrixDenseVectorProductBenchRelax<tags::Cell, float> CELLBMDVPBenchfloatRelax("CELL Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, float", 64ul*64ul*64ul, 100);
BandedMatrixDenseVectorProductBenchRelax<tags::Cell, double> CELLBMDVPBenchdoubleRelax("CELL Banded Matrix Dense Vector Product Relax Benchmark - matrix size: 64^3, double", 64ul*64ul*64ul, 100);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        DenseMatrixProductBench(const std::string & id, unsigned long size, unsigned long count) :
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
            for(unsigned long i = 0; i < _count; ++i)
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
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, float> MCSSEDMPBenchfloat2("MC SSE Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, double> MCSSEDMPBenchdouble2("MC SSE Matrix Product Benchmark dense/dense - matrix size: 256x256, double", 256, 10);
DenseMatrixProductBench<tags::CPU::SSE, float> SSEDMPBenchfloat3("SSE Matrix Product Benchmark dense/dense - matrix size: 2096x2096, float", 2096, 5);
DenseMatrixProductBench<tags::CPU::SSE, double> SSEDMPBenchdouble3("SSE Matrix Product Benchmark dense/dense - matrix size: 2096x2096, double", 2096, 5);
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, float> MCSSEDMPBenchfloat3("MC SSE Matrix Product Benchmark dense/dense - matrix size: 2096x2096, float", 2096, 5);
DenseMatrixProductBench<tags::CPU::MultiCore::SSE, double> MCSSEDMPBenchdouble3("MC SSE Matrix Product Benchmark dense/dense - matrix size: 2096x2096, double", 2096, 5);
#endif
#ifdef HONEI_CELL
DenseMatrixProductBench<tags::Cell, float> CELLDMPBenchfloat("CELL Matrix Product Benchmark dense/dense - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<tags::Cell, float> CELLDMPBenchfloat2("CELL Matrix Product Benchmark dense/dense - matrix size: 512x512, float", 512, 10);
DenseMatrixProductBench<tags::Cell, float> CELLDMPBenchfloat3("CELL Matrix Product Benchmark dense/dense - matrix size: 2048x2048, float", 2048, 5);
#endif

template <typename Tag_, typename DataType_>
class DenseMatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        DenseMatrixDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseMatrix<DataType_> dm0(_size, _size, DataType_(rand()));
            DenseVector<DataType_> dv0(_size, DataType_(rand()));
            for(unsigned long i = 0; i < _count; ++i)
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
#ifdef HONEI_CELL
DenseMatrixDenseVectorProductBench<tags::Cell, float> DMDVPBenchfloatCell("Cell: Matrix-Vector Product Benchmark dense/dense - matrix size: 64^2, float", 64ul*64, 10);
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        SparseMatrixProductBench(const std::string & id, unsigned long size, unsigned long count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            SparseMatrix<DataType_> sm(_size, _size, (unsigned long)((_size*_size)/10));
            for (typename SparseMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(rand());
                }
            }
            DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
            for(unsigned long i = 0; i < _count; ++i)
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
        unsigned long _count;
    public:
        BandedMatrixProductBench(const std::string & id, unsigned long size, unsigned long count) :
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
            bm.insert_band(1, dv);
            bm.insert_band(-1, dv);
            bm.insert_band(2, dv);
            bm.insert_band(-2, dv);
            bm.insert_band(5, dv);
            bm.insert_band(-5, dv);
            DenseMatrix<DataType_> dm(_size, _size, DataType_(rand()));
            for(unsigned long i = 0; i < _count; ++i)
            { 
                BENCHMARK(Product<Tag_>::value(bm, dm));
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm, dm));
            evaluate(info);
        }
};
BandedMatrixProductBench<tags::CPU, float> BMPBenchfloat2("Matrix Product Benchmark banded/dense - matrix size: 256x256, float", 256, 10);
BandedMatrixProductBench<tags::CPU, double> BMPBenchdouble2("Matrix Product Benchmark banded/dense - matrix size: 256x256, double", 256, 10);

template <typename Tag_, typename DataType_>
class SMELLDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
        std::string _file;
    public:
        SMELLDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count, std::string file) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
            _file = file;
        }

        virtual void run()
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _file;

            unsigned long non_zeros(MatrixIO<io_formats::M>::get_non_zeros(filename));
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DataType_> data(non_zeros);
            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            unsigned long rows, columns, ax, bx;
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrix<DataType_> tsmatrix(rows, columns, r, c, data);
            SparseMatrixELL<DataType_> smatrix(tsmatrix);

            DenseVector<DataType_> x(smatrix.rows());
            DenseVector<DataType_> y(smatrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DataType_(i) / 1.234;
            }

            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 1000 ; ++j)
                        {
                            Product<Tag_>::value(y, smatrix, x);
                        }
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            BenchmarkInfo info;
            info.flops = non_zeros * 2;
            evaluate(info * 1000);
        }
};
/// \todo Embed real world matrix with proper size
#ifdef HONEI_SSE
SMELLDenseVectorProductBench<tags::CPU::SSE, float> sse_SMELLDVPBenchfloat0("SSE SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10, "l8/area51_full_2.m");
SMELLDenseVectorProductBench<tags::CPU::SSE, double> sse_SMELLDVPBenchdouble0("SSE SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10, "l8/area51_full_2.m");
SMELLDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> mc_sse_SMELLDVPBenchfloat0("MC SSE SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10, "l8/area51_full_2.m");
SMELLDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> mc_sse_SMELLDVPBenchdouble0("MC SSE SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10, "l8/area51_full_2.m");
#endif
#ifdef HONEI_CUDA
SMELLDenseVectorProductBench<tags::GPU::CUDA, float> cudaSMELLDVPBenchfloat0("CUDA SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10, "l8/area51_full_2.m");
#ifdef HONEI_CUDA_DOUBLE
SMELLDenseVectorProductBench<tags::GPU::CUDA, double> cudaSMELLDVPBenchdouble0("CUDA SM 0 ELL Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10, "l8/area51_full_2.m");
#endif
#endif


template <typename Tag_, typename DataType_>
class Q1MatrixELLDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        Q1MatrixELLDenseVectorProductBench(const std::string & id, unsigned long size, unsigned long count) :
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

            bm1.insert_band(- (unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band(- (unsigned long)sqrt(_size) + 1, dv4.copy());
            bm1.insert_band(-1, dv4.copy());
            bm1.insert_band(0, dv4.copy());
            bm1.insert_band(1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size) - 1, dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size), dv4.copy());
            bm1.insert_band((unsigned long)sqrt(_size)+ 1, dv4.copy());

            //BandedMatrixQ1<DataType_> qm1(bm1);
            DenseVector<DataType_> dv2(_size, DataType_(4));
            DenseVector<DataType_> dv3(_size, DataType_(4));

            SparseMatrix<DataType_> sm(_size, _size, 10);
            for (unsigned long rows(0), cols(0) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(0), cols(1) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(1), cols(0) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(1024), cols(0) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(1025), cols(0) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(1026), cols(0) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(0), cols(1024) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(0), cols(1025) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }
            for (unsigned long rows(0), cols(1026) ; rows < _size && cols < _size; ++rows, ++cols)
            {
                sm(rows, cols) = DataType_(2);;
            }

            SparseMatrixELL<DataType_> smell(sm);
            for (unsigned long i(0) ; i < _count ; i++)
            {
                BENCHMARK(
                        for (unsigned long j(0) ; j < 10 ; ++j)
                        {
                            Product<Tag_>::value(dv3, smell, dv2);
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};
Q1MatrixELLDenseVectorProductBench<tags::CPU::SSE, float> sseQ1ELLDVPBenchfloat("SSE ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10);
Q1MatrixELLDenseVectorProductBench<tags::CPU::SSE, double> sseQ1ELLDVPBenchdouble("SSE ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10);
Q1MatrixELLDenseVectorProductBench<tags::CPU::MultiCore::SSE, float> mcsseQ1ELLDVPBenchfloat("MC SSE ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, float", 1025ul*1025, 10);
Q1MatrixELLDenseVectorProductBench<tags::CPU::MultiCore::SSE, double> mcsseQ1ELLDVPBenchdouble("MC SSE ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: L10, double", 1025ul*1025, 10);
#ifdef HONEI_CUDA
Q1MatrixELLDenseVectorProductBench<tags::GPU::CUDA, float> CUDAQ1ELLDVPBenchfloat("CUDA ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025, float",1025ul * 1025 , 10);
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixELLDenseVectorProductBench<tags::GPU::CUDA, double> CUDAQ1ELLDVPBenchdouble("CUDA ELL Matrix (Q1) Dense Vector Product Benchmark - matrix size: 1025*1025, double",1025ul * 1025 , 10);
#endif
#endif


template <typename DT_, typename Tag_>
class DenseMatrixDenseVectorProductSPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseMatrixDenseVectorProductSPUPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            int temp(Configuration::instance()->get_value("cell::product_dense_matrix_dense_vector_float", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < 60 ; k+=4)
                {
                    Configuration::instance()->set_value("cell::product_dense_matrix_dense_vector_float", j);
                    cores.push_back(stringify(j) +"SPUs" );
                    DenseVector<DT_> dv((k + 1) * 64, DT_(rand()));
                    DenseMatrix<DT_> dm((k + 1) * 64, (k + 1) * 64, DT_(rand()));

                    for(unsigned long i(0) ; i < 20 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                Product<Tag_>::value(dm, dv);
                                }
                                );
                    }
                    info = Product<>::get_benchmark_info(dm, dv);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::product_dense_matrix_dense_vector_float", temp);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifdef HONEI_CELL
DenseMatrixDenseVectorProductSPUPlot<float, tags::Cell> DMDVPSPUF("Cell Dense Matrix Dense Vector Product Benchmark - SPU Count: 1 to 6 - float");
#endif

template <typename DT_, typename Tag_>
class BandedMatrixProductQ1SPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        BandedMatrixProductQ1SPUPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            int temp(Configuration::instance()->get_value("cell::product_banded_matrix_dense_vector_float", 4));
            int temp2(Configuration::instance()->get_value("cell::product_banded_matrix_dense_vector_double", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < 152 ; k+=5)
                {
                    Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_float", j);
                    Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_double", j);
                    cores.push_back(stringify(j) +"SPUs" );
                    unsigned long _size((k+1) * 8192);
                    DenseVector<DT_> dv0((k + 1) * 8192, DT_(rand()));
                    DenseVector<DT_> dv1((k + 1) * 8192, DT_(rand()));
                    BandedMatrix<DT_> bm((k + 1) * 8192, dv1);
                    bm.insert_band(- (unsigned long)sqrt(_size) - 1, dv1.copy());
                    bm.insert_band(- (unsigned long)sqrt(_size), dv1.copy());
                    bm.insert_band(- (unsigned long)sqrt(_size) + 1, dv1.copy());
                    bm.insert_band(-1, dv1.copy());
                    bm.insert_band(0, dv1.copy());
                    bm.insert_band(1, dv1.copy());
                    bm.insert_band((unsigned long)sqrt(_size) - 1, dv1.copy());
                    bm.insert_band((unsigned long)sqrt(_size), dv1.copy());
                    bm.insert_band((unsigned long)sqrt(_size)+ 1, dv1.copy());

                    for(unsigned long i(0) ; i < 10 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                Product<Tag_>::value(bm, dv0);
                                }
                                );
                    }
                    info = Product<>::get_benchmark_info(bm, dv0);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_float", temp);
            Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_double", temp2);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 10);
        }
};
#ifdef HONEI_CELL
BandedMatrixProductQ1SPUPlot<float, tags::Cell> BMPQ1SPUF("Cell Q1 Banded Matrix Dense Vector Product Benchmark - SPU Count: 1 to 6 - float");
BandedMatrixProductQ1SPUPlot<double, tags::Cell> BMPQ1SPUD("Cell Q1 Banded Matrix Dense Vector Product Benchmark - SPU Count: 1 to 6 - double");
#endif


template <typename DT_, typename Tag_>
class BandedMatrixProductRelaxSPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        BandedMatrixProductRelaxSPUPlot(const std::string & id, int x) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _plots = true;
            _x = x;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            int temp(Configuration::instance()->get_value("cell::product_banded_matrix_dense_vector_float", 4));
            int temp2(Configuration::instance()->get_value("cell::product_banded_matrix_dense_vector_double", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < _x ; k+=20)
                {
                    Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_float", j);
                    Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_double", j);
                    cores.push_back(stringify(j) +"SPUs" );

                    DenseVector<DT_> dv0((k + 1) * 8192, DT_(rand()));
                    DenseVector<DT_> dv1((k + 1) * 8192, DT_(rand()));
                    BandedMatrix<DT_> bm((k + 1) * 8192, dv1);
                    bm.insert_band(3, dv1.copy());
                    bm.insert_band(-3, dv1.copy());
                    bm.insert_band(15, dv1.copy());

                    for(unsigned long i(0) ; i < 10 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                Product<Tag_>::value(bm, dv0);
                                }
                                );
                    }
                    info = Product<>::get_benchmark_info(bm, dv0);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_float", temp);
            Configuration::instance()->set_value("cell::product_banded_matrix_dense_vector_double", temp2);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 10);
        }
};
#ifdef HONEI_CELL
BandedMatrixProductRelaxSPUPlot<float, tags::Cell> BMPRelaxSPUF("Cell Relax Banded Matrix Dense Vector Product Benchmark - SPU Count: 1 to 6 - float", 400);
BandedMatrixProductRelaxSPUPlot<double, tags::Cell> BMPRelaxSPUD("Cell Relax Banded Matrix Dense Vector Product Benchmark - SPU Count: 1 to 6 - double", 180);
#endif

template <typename DT_>
class BandedMatrixProductQ1VSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        BandedMatrixProductQ1VSPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(tags::CPU::SSE::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            // mc::sse
            for (unsigned long j(1) ; j < 160 ; j+=10)
            {
                cores.push_back(tags::CPU::MultiCore::SSE::name);
                unsigned long _size((j+1) * 8192);
                DenseVector<DT_> dv0((j + 1) * 8192, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 8192, DT_(rand()));
                BandedMatrix<DT_> bm((j + 1) * 8192, dv1);
                bm.insert_band(- (unsigned long)sqrt(_size) - 1, dv1.copy());
                bm.insert_band(- (unsigned long)sqrt(_size), dv1.copy());
                bm.insert_band(- (unsigned long)sqrt(_size) + 1, dv1.copy());
                bm.insert_band(-1, dv1.copy());
                bm.insert_band(0, dv1.copy());
                bm.insert_band(1, dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size) - 1, dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size), dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size)+ 1, dv1.copy());

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::MultiCore::SSE>::value(bm, dv0));
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long j(1) ; j < 160 ; j+=10)
            {
                cores.push_back(tags::CPU::SSE::name);
                unsigned long _size((j+1) * 8192);
                DenseVector<DT_> dv0((j + 1) * 8192, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 8192, DT_(rand()));
                BandedMatrix<DT_> bm((j + 1) * 8192, dv1);
                bm.insert_band(- (unsigned long)sqrt(_size) - 1, dv1.copy());
                bm.insert_band(- (unsigned long)sqrt(_size), dv1.copy());
                bm.insert_band(- (unsigned long)sqrt(_size) + 1, dv1.copy());
                bm.insert_band(-1, dv1.copy());
                bm.insert_band(0, dv1.copy());
                bm.insert_band(1, dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size) - 1, dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size), dv1.copy());
                bm.insert_band((unsigned long)sqrt(_size)+ 1, dv1.copy());

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::SSE>::value(bm, dv0));
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
BandedMatrixProductQ1VSPlot<float> BMPQ1VSF("MC vs SSE Q1-BandedMatrix Product Benchmark - float");
BandedMatrixProductQ1VSPlot<double> BMPQ1VSD("MC vs SSE Q1-BandedMatrix Product Benchmark - double");
#endif

template <typename DT_>
class BandedMatrixProductRelaxVSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        BandedMatrixProductRelaxVSPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(tags::CPU::SSE::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            // mc::sse
            for (unsigned long j(1) ; j < 160 ; j+=10)
            {
                cores.push_back(tags::CPU::MultiCore::SSE::name);
                unsigned long _size((j + 1) * 8192);
                DenseVector<DT_> dv0(_size, DT_(rand()));
                DenseVector<DT_> dv1(_size, DT_(rand()));
                BandedMatrix<DT_> bm(_size, dv1);
                bm.insert_band(3, dv1.copy());
                bm.insert_band(-3, dv1.copy());
                bm.insert_band(15, dv1.copy());

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::MultiCore::SSE>::value(bm, dv0));
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long j(1) ; j < 160 ; j+=10)
            {
                cores.push_back(tags::CPU::SSE::name);
                unsigned long _size((j+1) * 8192);
                DenseVector<DT_> dv0(_size, DT_(rand()));
                DenseVector<DT_> dv1(_size, DT_(rand()));
                BandedMatrix<DT_> bm(_size, dv1);
                bm.insert_band(3, dv1.copy());
                bm.insert_band(-3, dv1.copy());
                bm.insert_band(15, dv1.copy());

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::SSE>::value(bm, dv0));
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
BandedMatrixProductRelaxVSPlot<float> BMPRVSF("MC vs SSE Relax-BandedMatrix Product Benchmark - float");
BandedMatrixProductRelaxVSPlot<double> BMPRVSD("MC vs SSE Relax-BandedMatrix Product Benchmark - double");
#endif

template <typename DT_>
class DenseMatrixProductVSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseMatrixProductVSPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(tags::CPU::SSE::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            // mc::sse
            for (unsigned long k(1) ; k < 69 ; k+=5)
            {
                cores.push_back(tags::CPU::MultiCore::SSE::name);
                DenseMatrix<DT_> dm0((k + 1) * 40, (k + 1) * 40, DT_(rand()));
                DenseMatrix<DT_> dm1((k + 1) * 40, (k + 1) * 40, DT_(rand()));

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::MultiCore::SSE>::value(dm0, dm1));
                }
                info = Product<>::get_benchmark_info(dm0, dm1);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long k(1) ; k < 69 ; k+=5)
            {
                cores.push_back(tags::CPU::SSE::name);
                DenseMatrix<DT_> dm0((k + 1) * 40, (k + 1) * 40, DT_(rand()));
                DenseMatrix<DT_> dm1((k + 1) * 40, (k + 1) * 40, DT_(rand()));

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::SSE>::value(dm0, dm1));
                }
                info = Product<>::get_benchmark_info(dm0, dm1);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
DenseMatrixProductVSPlot<float> DMPVSF("MC vs SSE DenseMatrix Product Benchmark - float");
DenseMatrixProductVSPlot<double> DMPVSD("MC vs SSE DenseMatrix Product Benchmark - double");
#endif


template <typename DT_>
class DenseMatrixDenseVectorProductVSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseMatrixDenseVectorProductVSPlot(const std::string & id) :
            Benchmark(id)
        {
            register_tag(tags::CPU::SSE::name);
            _plots = true;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;

            // mc::sse
            for (unsigned long j(1) ; j < 95 ; j+=5)
            {
                cores.push_back(tags::CPU::MultiCore::SSE::name);
                DenseVector<DT_> dv((j + 1) * 64, DT_(rand()));
                DenseMatrix<DT_> dm((j + 1) * 64, (j + 1) * 64, DT_(rand()));

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::MultiCore::SSE>::value(dm, dv));
                }
                info = Product<>::get_benchmark_info(dm, dv);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long j(1) ; j < 95 ; j+=5)
            {
                cores.push_back(tags::CPU::SSE::name);
                DenseVector<DT_> dv((j + 1) * 64, DT_(rand()));
                DenseMatrix<DT_> dm((j + 1) * 64, (j + 1) * 64, DT_(rand()));

                for(unsigned long i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(Product<tags::CPU::SSE>::value(dm, dv));
                }
                info = Product<>::get_benchmark_info(dm, dv);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
DenseMatrixDenseVectorProductVSPlot<float> DMDVPVSF("MC vs SSE DenseMatrix DenseVector Product Benchmark - float");
DenseMatrixDenseVectorProductVSPlot<double> DMDVPVSD("MC vs SSE DenseMatrix DenseVector Product Benchmark - double");
#endif
