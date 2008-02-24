/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/libla/reduction.hh>

using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class DenseMatrixRowSumVectorBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DenseMatrixRowSumVectorBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size);
            DenseMatrix<DataType_> dm(_size, _size, static_cast<DataType_>(rand()));
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(dv = (Reduction<rt_sum, Tag_>::value(dm)));
            }
            BenchmarkInfo info(Reduction<rt_sum>::get_benchmark_info(dm));
            evaluate(info);
        }
};
DenseMatrixRowSumVectorBench<tags::CPU, float> DMRSVBenchfloat1("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU, double> DMRSVBenchdouble1("Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore, float> DMRSVBenchfloat2("MC: Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore, double> DMRSVBenchdouble2("MC: Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
DenseMatrixRowSumVectorBench<tags::CPU::SSE, float> SSEDMRSVBenchfloat1("SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::SSE, double> SSEDMRSVBenchdouble1("SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, float> SSEDMRSVBenchfloat2("MC::SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, double> SSEDMRSVBenchdouble2("MC::SSE Dense Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#endif

template <typename Tag_, typename DataType_>
class SparseMatrixRowSumVectorBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        SparseMatrixRowSumVectorBench(const std::string & id, unsigned long size, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            DenseVector<DataType_> dv(_size);
            SparseMatrix<DataType_> sm(_size, _size, (unsigned long)(_size/10));
            for (typename MutableMatrix<DataType_>::ElementIterator i_end(sm.end_elements()), i(sm.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DataType_(rand());
                }
            }
            for(int i = 0; i < _count; ++i)
            {
                BENCHMARK(dv = (Reduction<rt_sum, Tag_>::value(sm)));
            }
            BenchmarkInfo info(Reduction<rt_sum>::get_benchmark_info(sm));
            evaluate(info);
        }
};
SparseMatrixRowSumVectorBench<tags::CPU, float> SMRSVBenchfloat1("Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU, double> SMRSVBenchdouble1("Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore, float> SMRSVBenchfloat2("MC: Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore, double> SMRSVBenchdouble2("MC: Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#ifdef HONEI_SSE
SparseMatrixRowSumVectorBench<tags::CPU::SSE, float> SSESMRSVBenchfloat1("SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::SSE, double> SSESMRSVBenchdouble1("SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, float> SSESMRSVBenchfloat2("MC::SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, float", 4096, 10);
SparseMatrixRowSumVectorBench<tags::CPU::MultiCore::SSE, double> SSESMRSVBenchdouble2("MC::SSE Sparse Matrix Row Sum Vector Benchmark - matrix size: 4096x4096, double", 4096, 10);
#endif

template <typename DT_, typename Tag_>
class DenseVectorReductionSPUPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseVectorReductionSPUPlot(const std::string & id) :
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

            int temp(Configuration::instance()->get_value("cell::reduction_sum_dense_float", 4));
            int temp2(Configuration::instance()->get_value("cell::reduction_sum_dense_double", 4));

            int max_spu(6);

            for (unsigned long j(1) ; j <= max_spu ; ++j)
            {
                for (unsigned long k(1) ; k < 81 ; k+=5)
                {
                    Configuration::instance()->set_value("cell::reduction_sum_dense_float", j);
                    Configuration::instance()->set_value("cell::reduction_sum_dense_double", j);
                    cores.push_back(stringify(j) +"SPUs" );
                    DenseVector<DT_> dv0(k * 150000, DT_(rand()));
                    DenseVector<DT_> dv1(k * 150000, DT_(rand()));

                    for(int i(0) ; i < 20 ; ++i)
                    {
                        BENCHMARK(
                                for (unsigned long l(0) ; l < 5 ; ++l)
                                {
                                (Reduction<rt_sum, Tag_>::value(dv0));
                                }
                                );
                    }
                    info = Reduction<rt_sum>::get_benchmark_info(dv0);
                    infolist.push_back(info * 5);
                    std::cout << ".";
                    std::cout.flush();
                }
            }
            Configuration::instance()->set_value("cell::reduction_sum_dense_float", temp);
            Configuration::instance()->set_value("cell::reduction_sum_dense_double", temp2);
            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifdef HONEI_CELL
DenseVectorReductionSPUPlot<float, tags::Cell> DVRSPUF("Cell Dense Vector Reduction to Sum Benchmark - SPU Count: 1 to 6 - float");
DenseVectorReductionSPUPlot<double, tags::Cell> DVRSPUD("Cell Dense Vector Recution to Sum Benchmark - SPU Count: 1 to 6 - double");
#endif


template <typename DT_>
class DenseVectorReductionVSPlot :
    public Benchmark
{
    private:
        int _x;

    public:
        DenseVectorReductionVSPlot(const std::string & id) :
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
                DenseVector<DT_> dv((j + 1) * 131072, DT_(rand()));

                for(int i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK((Reduction<rt_sum, tags::CPU::MultiCore::SSE>::value(dv)));
                }
                info = Reduction<rt_sum>::get_benchmark_info(dv);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            // sse
            for (unsigned long j(1) ; j < 95 ; j+=5)
            {
                cores.push_back(tags::CPU::SSE::name);
                DenseVector<DT_> dv((j + 1) * 131072, DT_(rand()));

                for(int i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK((Reduction<rt_sum, tags::CPU::SSE>::value(dv)));
                }
                info = Reduction<rt_sum>::get_benchmark_info(dv);
                infolist.push_back(info);
                std::cout<<".";
                std::cout.flush();
            }

            std::cout<<std::endl;
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifdef HONEI_SSE
DenseVectorReductionVSPlot<float> DVRVSF("MC vs SSE DenseVector Reduction to sum Benchmark - float");
DenseVectorReductionVSPlot<double> DVRVSD("MC vs SSE DenseVector Reduction to sum Benchmark - double");
#endif
