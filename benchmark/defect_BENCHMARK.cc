/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

#include <string>
#endif

#include <honei/math/defect.hh>

using namespace honei;

template <typename Tag_, typename DataType_>
class DefectBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
    public:
        DefectBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> rhs(_size, DataType_(2.3456789));
            DenseVector<DataType_> x(_size, DataType_(0.123456789));

            DenseVector<DataType_> omni_band(_size, DataType_(4711));
            BandedMatrixQx<Q1Type, DataType_> A(_size, omni_band.copy(), omni_band.copy(),omni_band.copy(),omni_band.copy(),omni_band.copy(),omni_band.copy(),omni_band.copy(),omni_band.copy(),omni_band.copy());

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        for(unsigned long j(0) ; j < 100 ; ++j)
                        {
                        Defect<Tag_>::value(rhs, A, x);
                        });
            }

            evaluate();
        }
};

DefectBench<tags::CPU, float> defect_bench_float("Defect bench float", 16641, 100);
DefectBench<tags::CPU, double> defect_bench_double("Defect bench double", 16641, 100);
#ifdef HONEI_SSE
DefectBench<tags::CPU::SSE, float> defect_bench_float_sse("Defect bench float SSE", 16641, 100);
DefectBench<tags::CPU::SSE, double> defect_bench_double_sse("Defect bench double SSE", 16641, 100);
#endif
#ifdef HONEI_CELL
//DefectBench<tags::Cell, float> defect_bench_float_cell("Defect bench float Cell", 16641, 100);
//DefectBench<tags::Cell, double> defect_bench_double_cell("Defect bench double Cell", 16641, 100);
#endif
#ifdef HONEI_CUDA
DefectBench<tags::GPU::CUDA, float> defect_bench_float_cuda("Defect bench float CUDA", 16641, 100);
#endif
