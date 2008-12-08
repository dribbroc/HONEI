/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/la/product.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>
class Q1MatrixDenseVectorProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        Q1MatrixDenseVectorProductBench(const std::string & id, unsigned long size, int count) :
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
                            cuda_thread_synchronize();
#endif
                        }
                        );
            }
            BenchmarkInfo info(Product<>::get_benchmark_info(bm1, dv2));
            evaluate(info * 10);
        }
};

Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat1("CPU: Size: 1089, float", 1089ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat2("CPU: Size: 4225, float", 4225ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat3("CPU: Size: 16641, float", 16641ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat4("CPU: Size: 66049, float", 66049ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat5("CPU: Size: 263169, float", 263169ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat6("CPU: Size: 1050625, float", 1050625ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU, float> Q1DVPBenchfloat7("CPU: Size: 2198401, float", 2198401ul, 10);

#ifdef HONEI_SSE
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat1("SSE: Size: 1089, float", 1089ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat2("SSE: Size: 4225, float", 4225ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat3("SSE: Size: 16641, float", 16641ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat4("SSE: Size: 66049, float", 66049ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat5("SSE: Size: 263169, float", 263169ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat6("SSE: Size: 1050625, float", 1050625ul, 10);
Q1MatrixDenseVectorProductBench<tags::CPU::SSE, float> sseQ1DVPBenchfloat7("SSE: Size: 2198401, float", 2198401ul, 10);
#endif

#ifdef HONEI_CELL
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat1("Cell: Size: 1089, float", 1089ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat2("Cell: Size: 4225, float", 4225ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat3("Cell: Size: 16641, float", 16641ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat4("Cell: Size: 66049, float", 66049ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat5("Cell: Size: 263169, float", 263169ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat6("Cell: Size: 1050625, float", 1050625ul, 10);
Q1MatrixDenseVectorProductBench<tags::Cell, float> cellQ1DVPBenchfloat7("Cell: Size: 2198401, float", 2198401ul, 10);
#endif

#ifdef HONEI_CUDA
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat1("CUDA: Size: 1089, float", 1089ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat2("CUDA: Size: 4225, float", 4225ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat3("CUDA: Size: 16641, float", 16641ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat4("CUDA: Size: 66049, float", 66049ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat5("CUDA: Size: 263169, float", 263169ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat6("CUDA: Size: 1050625, float", 1050625ul, 10);
Q1MatrixDenseVectorProductBench<tags::GPU::CUDA, float> cudaQ1DVPBenchfloat7("CUDA: Size: 2198401, float", 2198401ul, 10);
#endif
