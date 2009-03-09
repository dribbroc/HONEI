/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/scaled_product_sum_norm.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>


using namespace std;
using namespace honei;


template <typename Tag_, typename DataType_>
class SPSNBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SPSNBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
                DenseVector<DataType_> x(_size, static_cast<DataType_>(rand()));
                DenseVector<DataType_> y(_size, static_cast<DataType_>(rand()));

                DenseVector<DataType_> ll(x.copy());
                DenseVector<DataType_> ld(x.copy());
                DenseVector<DataType_> lu(x.copy());
                DenseVector<DataType_> dl(x.copy());
                DenseVector<DataType_> dd(x.copy());
                DenseVector<DataType_> du(x.copy());
                DenseVector<DataType_> ul(x.copy());
                DenseVector<DataType_> ud(x.copy());
                DenseVector<DataType_> uu(x.copy());

                BandedMatrixQ1<DataType_> A(_size, ll, ld, lu, dl, dd, du, ul, ud, uu);

                DataType_ a(.123456789);
                DataType_ b(1234.56789);

                DataType_ result(0.);
                for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long l(0) ; l < 10 ; ++l)
                        {
                        result = ScaledProductSumNorm<Tag_>::value(a, y, b, A, x);
#ifdef HONEI_CUDA
                        //cuda_thread_synchronize();
#endif
                        }
                        );
            }
            BenchmarkInfo info(ScaledProductSumNorm<tags::CPU>::get_benchmark_info(x));
            evaluate(info * 10);
        }
};
#ifdef HONEI_SSE
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m6("SSE SPSN   Benchmark - vector size: 27, float", 27ul, 10000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m5("SSE SPSN   Benchmark - vector size: 81, float", 81ul, 10000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m4("SSE SPSN   Benchmark - vector size: 289, float", 289ul, 1000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m3("SSE SPSN   Benchmark - vector size: 1089, float", 1089ul, 1000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m2("SSE SPSN   Benchmark - vector size: 4225, float", 4225ul, 1000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_m1("SSE SPSN   Benchmark - vector size: 16641, float", 16641ul, 1000);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_0("SSE SPSN   Benchmark - vector size: 66049, float", 66049ul, 100);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_1("SSE SPSN   Benchmark - vector size: 263169, float", 263169ul, 100);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_2("SSE SPSN   Benchmark - vector size: 1050625, float", 1050625ul, 100);
SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_3("SSE SPSN   Benchmark - vector size: 4198401, float", 4198401ul, 100);
//SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_4("SSE SPSN   Benchmark - vector size: 16785409, float", 16785409ul, 10);
//SPSNBench<tags::CPU::SSE, float> spsn_opt_float_sse_5("SSE SPSN   Benchmark - vector size: 25130169, float", 25130169ul, 10);
#endif

template <typename Tag_, typename DataType_>
class SPSNBench_TUTORIAL :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SPSNBench_TUTORIAL(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
                DenseVector<DataType_> x(_size, static_cast<DataType_>(rand()));
                DenseVector<DataType_> y(_size, static_cast<DataType_>(rand()));

                DenseVector<DataType_> ll(x.copy());
                DenseVector<DataType_> ld(x.copy());
                DenseVector<DataType_> lu(x.copy());
                DenseVector<DataType_> dl(x.copy());
                DenseVector<DataType_> dd(x.copy());
                DenseVector<DataType_> du(x.copy());
                DenseVector<DataType_> ul(x.copy());
                DenseVector<DataType_> ud(x.copy());
                DenseVector<DataType_> uu(x.copy());

                BandedMatrixQ1<DataType_> A(_size, ll, ld, lu, dl, dd, du, ul, ud, uu);

                DataType_ a(.123456789);
                DataType_ b(1234.56789);

                DataType_ result(0.);
                for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(
                        for (unsigned long l(0) ; l < 10 ; ++l)
                        {
                        result = ScaledProductSumNorm_TUTORIAL<Tag_>::value(a, y, b, A, x);
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
                        }
                        );
            }
            BenchmarkInfo info(ScaledProductSumNorm_TUTORIAL<tags::CPU>::get_benchmark_info(x));
            evaluate(info * 10);
        }
};
#ifdef HONEI_CUDA
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m6("CUDA SPSN   Benchmark - vector size: 27, float", 27ul, 100);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m5("CUDA SPSN   Benchmark - vector size: 81, float", 81ul, 100);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m4("CUDA SPSN   Benchmark - vector size: 289, float", 289ul, 1000);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m3("CUDA SPSN   Benchmark - vector size: 1089, float", 1089ul, 1000);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m2("CUDA SPSN   Benchmark - vector size: 4225, float", 4225ul, 1000);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_m1("CUDA SPSN   Benchmark - vector size: 16641, float", 16641ul, 1000);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_0("CUDA SPSN Benchmark - vector size: 66049, float", 66049ul, 100);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_1("CUDA SPSN Benchmark - vector size: 263169, float", 263169ul, 100);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_2("CUDA SPSN Benchmark - vector size: 1050625, float", 1050625ul, 100);
SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_3("CUDA SPSN Benchmark - vector size:  4198401, float", 4198401ul, 100);
//SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_4("CUDA SPSN Benchmark - vector size: 16785409, float", 16785409ul, 10);
//SPSNBench_TUTORIAL<tags::GPU::CUDA, float> spsn_tut_float_cuda_5("CUDA SPSN Benchmark - vector size: 25130169, float", 25130169ul, 10);
#endif
