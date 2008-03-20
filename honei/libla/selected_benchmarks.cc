/* vim: set sw=4 sts=4 et nofoldenable : */
#include <benchmark/benchmark.cc>

#include <honei/libla/sum.hh>
#include <honei/libla/product.hh>
#include <honei/libla/reduction.hh>
#include <honei/libla/dot_product.hh>
#include <honei/libla/scaled_sum.hh>

using namespace honei;


template <typename DT_, typename Tag_>
class DenseVectorRTSBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DenseVectorRTSBench(const std::string & id, int count) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _plots = true;
            _count = count;
        }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv((j + 1) * 131072, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            (Reduction<rt_sum, Tag_>::value(dv));
                            }
                            );
                }
                info = Reduction<rt_sum>::get_benchmark_info(dv);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//DenseVectorRTSBench<float, tags::CPU> DVRBF("SingleCore DenseVector Reduction to Sum Benchmark - float", 20);
//DenseVectorRTSBench<double, tags::CPU> DVRBD("SingleCore DenseVector Reduction to Sum Benchmark - double", 20);
//DenseVectorRTSBench<float, tags::CPU::MultiCore> DVRBMCF("MultiCore DenseVector Reduction to Sum Benchmark - float", 20);
//DenseVectorRTSBench<double, tags::CPU::MultiCore> DVRBMCD("MultiCore DenseVector Reduction to Sum Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorRTSBench<float, tags::CPU::SSE> DVRBSSEF("SSE DenseVector Reduction to Sum Benchmark - float", 90);
DenseVectorRTSBench<double, tags::CPU::SSE> DVRBSSED("SSE DenseVector Reduction to Sum Benchmark - double", 90);
DenseVectorRTSBench<float, tags::CPU::MultiCore::SSE> DVRBMCSSEF("MultiCore SSE DenseVector Reduction to Sum Benchmark - float", 90);
DenseVectorRTSBench<double, tags::CPU::MultiCore::SSE> DVRBMCSSED("MultiCore SSE DenseVector Reduction to Sum Benchmark - double", 90);
#endif
#elif HONEI_CELL
DenseVectorRTSBench<float, tags::Cell> DVRBCF("CELL DenseVector Reduction to Sum Benchmark - float", 90);
DenseVectorRTSBench<double, tags::Cell> DVRBCD("CELL DenseVector Reduction to Sum Benchmark - double", 90);
#endif


template <typename DT_, typename Tag_>
class DenseVectorDotProductBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DenseVectorDotProductBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            DotProduct<Tag_>::value(dv0, dv1);
                            }
                            );
                }
                info = DotProduct<>::get_benchmark_info(dv0, dv1);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//DenseVectorDotProductBench<float, tags::CPU> DVDPBF("SingleCore DenseVector DotProduct Benchmark - float", 20);
//DenseVectorDotProductBench<double, tags::CPU> DVDPBD("SingleCore DenseVector DotProduct Benchmark - double", 20);
//DenseVectorDotProductBench<float, tags::CPU::MultiCore> DVDPBMCF("MultiCore DenseVector DotProduct Benchmark - float", 20);
//DenseVectorDotProductBench<double, tags::CPU::MultiCore> DVDPBMCD("MultiCore DenseVector DotProduct Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorDotProductBench<float, tags::CPU::SSE> DVDPBSSEF("SSE DenseVector DotProduct Benchmark - float", 90);
DenseVectorDotProductBench<double, tags::CPU::SSE> DVDPBSSED("SSE DenseVector DotProduct Benchmark - double", 90);
DenseVectorDotProductBench<float, tags::CPU::MultiCore::SSE> DVDPBMCSSEF("MultiCore SSE DenseVector DotProduct Benchmark - float", 90);
DenseVectorDotProductBench<double, tags::CPU::MultiCore::SSE> DVDPBMCSSED("MultiCore SSE DenseVector DotProduct Benchmark - double", 90);
#endif
#elif HONEI_CELL
DenseVectorDotProductBench<float, tags::Cell> DVDPBCF("CELL DenseVector DotProduct Benchmark - float", 90);
DenseVectorDotProductBench<double, tags::Cell> DVDPBCD("CELL DenseVector DotProduct Benchmark - double", 90);
#endif


template <typename DT_, typename Tag_>
class DenseVectorSumBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DenseVectorSumBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            Sum<Tag_>::value(dv0, dv1);
                            }
                            );
                }
                info = Sum<>::get_benchmark_info(dv0, dv1);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//DenseVectorSumBench<float, tags::CPU> DVSBF("SingleCore DenseVector Sum Benchmark - float", 20);
//DenseVectorSumBench<double, tags::CPU> DVSBD("SingleCore DenseVector Sum Benchmark - double", 20);
//DenseVectorSumBench<float, tags::CPU::MultiCore> DVSBMCF("MultiCore DenseVector Sum Benchmark - float", 20);
//DenseVectorSumBench<double, tags::CPU::MultiCore> DVSBMCD("MultiCore DenseVector Sum Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorSumBench<float, tags::CPU::SSE> DVSBSSEF("SSE DenseVector Sum Benchmark - float", 90);
DenseVectorSumBench<double, tags::CPU::SSE> DVSBSSED("SSE DenseVector Sum Benchmark - double", 90);
DenseVectorSumBench<float, tags::CPU::MultiCore::SSE> DVSBMCSSEF("MultiCore SSE DenseVector Sum Benchmark - float", 90);
DenseVectorSumBench<double, tags::CPU::MultiCore::SSE> DVSBMCSSED("MultiCore SSE DenseVector Sum Benchmark - double", 90);
#endif
#elif HONEI_CELL
DenseVectorSumBench<float, tags::Cell> DVSBCF("CELL DenseVector Sum Benchmark - float", 90);
DenseVectorSumBench<double, tags::Cell> DVSBCD("CELL DenseVector Sum Benchmark - double", 90);
#endif


template <typename DT_, typename Tag_>
class DenseVectorScaledSumBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DenseVectorScaledSumBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                DT_ alpha(rand());
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            ScaledSum<Tag_>::value(dv0, dv1, alpha);
                            }
                            );
                }
                info = ScaledSum<>::get_benchmark_info(dv0, dv1, alpha);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//DenseVectorScaledSumBench<float, tags::CPU> DVSSBF("SingleCore DenseVector ScaledSum Benchmark - float", 20);
//DenseVectorScaledSumBench<double, tags::CPU> DVSSBD("SingleCore DenseVector ScaledSum Benchmark - double", 20);
//DenseVectorScaledSumBench<float, tags::CPU::MultiCore> DVSSBMCF("MultiCore DenseVector ScaledSum Benchmark - float", 20);
//DenseVectorScaledSumBench<double, tags::CPU::MultiCore> DVSSBMCD("MultiCore DenseVector ScaledSum Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorScaledSumBench<float, tags::CPU::SSE> DVSSBSSEF("SSE DenseVector ScaledSum Benchmark - float", 90);
DenseVectorScaledSumBench<double, tags::CPU::SSE> DVSSBSSED("SSE DenseVector ScaledSum Benchmark - double", 90);
DenseVectorScaledSumBench<float, tags::CPU::MultiCore::SSE> DVSSBMCSSEF("MultiCore SSE DenseVector ScaledSum Benchmark - float", 90);
DenseVectorScaledSumBench<double, tags::CPU::MultiCore::SSE> DVSSBMCSSED("MultiCore SSE DenseVector ScaledSum Benchmark - double", 90);
#endif
#elif HONEI_CELL
DenseVectorScaledSumBench<float, tags::Cell> DVSSBCF("CELL DenseVector ScaledSum Benchmark - float", 90);
DenseVectorScaledSumBench<double, tags::Cell> DVSSBCD("CELL DenseVector ScaledSum Benchmark - double", 90);
#endif

template <typename DT_, typename Tag_>
class DMDVProductBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DMDVProductBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv((j + 1) * 64, DT_(rand()));
                DenseMatrix<DT_> dm((j + 1) * 64, (j + 1) * 64, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            Product<Tag_>::value(dm, dv);
                            }
                            );
                }
                info = Product<>::get_benchmark_info(dm, dv);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//DMDVProductBench<float, tags::CPU> DMDVPBTPF("SingleCore DenseMatrix DenseVector Product Benchmark - float", 30);
//DMDVProductBench<double, tags::CPU> DMDVPBTPD("SingleCore DenseMatrix DenseVector Product Benchmark - double", 30);
//DMDVProductBench<float, tags::CPU::MultiCore> DMDVPBMCF("MultiCore DenseMatrix DenseVector Product Benchmark - float", 30);
//DMDVProductBench<double, tags::CPU::MultiCore> DMDVPBMCD("MultiCore DenseMatrix DenseVector Product Benchmark - double", 30);
#ifdef HONEI_SSE
DMDVProductBench<float, tags::CPU::SSE> DMDVPBSSEF("SSE DenseMatrix DenseVector Product Benchmark - float", 90);
DMDVProductBench<double, tags::CPU::SSE> DMDVPBSSED("SSE DenseMatrix DenseVector Product Benchmark - double", 90);
DMDVProductBench<float, tags::CPU::MultiCore::SSE> DMDVPBMCSSEF("MultiCore SSE DenseMatrix DenseVector Product Benchmark - float", 90);
DMDVProductBench<double, tags::CPU::MultiCore::SSE> DMDVPBMCSSED("MultiCore SSE DenseMatrix DenseVector Product Benchmark - double", 90);
#endif
#elif HONEI_CELL
DMDVProductBench<float, tags::Cell> DMDVPBCF("CELL DenseMatrix DenseVector Product Benchmark - float", 54);
DMDVProductBench<double, tags::Cell> DMDVPBCD("CELL DenseMatrix DenseVector Product Benchmark - double", 20);
#endif

template <typename DT_, typename Tag_>
class DenseMatrixProductBench :
    public Benchmark
{
    private:
        int _count;

    public:
        DenseMatrixProductBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
                DenseMatrix<DT_> dm0((j + 1) * 40, (j + 1) * 40, DT_(rand()));
                DenseMatrix<DT_> dm1((j + 1) * 40, (j + 1) * 40, DT_(rand()));
                for(int i(0) ; i < 5 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 5 ; ++k)
                            {
                            Product<Tag_>::value(dm0, dm1);
                            }
                            );
                }
                info = Product<>::get_benchmark_info(dm0, dm1);
                infolist.push_back(info * 5);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 5);
        }
};
#ifndef HONEI_CELL
//DenseMatrixProductBench<float, tags::CPU> DMPBTPF("SingleCore DenseMatrix Product Benchmark - float", 15);
//DenseMatrixProductBench<double, tags::CPU> DMPBTPD("SingleCore DenseMatrix Product Benchmark - double", 15);
//DenseMatrixProductBench<float, tags::CPU::MultiCore> DMPBMCF("MultiCore DenseMatrix Product Benchmark - float", 15);
//DenseMatrixProductBench<double, tags::CPU::MultiCore> DMPBMCD("MultiCore DenseMatrix Product Benchmark - double", 15);
#ifdef HONEI_SSE
DenseMatrixProductBench<float, tags::CPU::SSE> DMPBSSEF("SSE DenseMatrix Product Benchmark - float", 54);
DenseMatrixProductBench<double, tags::CPU::SSE> DMPBSSED("SSE DenseMatrix Product Benchmark - double", 54);
DenseMatrixProductBench<float, tags::CPU::MultiCore::SSE> DMPBMCSSEF("MultiCore SSE DenseMatrix Product Benchmark - float", 54);
DenseMatrixProductBench<double, tags::CPU::MultiCore::SSE> DMPBMCSSED("MultiCore SSE DenseMatrix Product Benchmark - double", 54);
#endif
#elif HONEI_CELL
DenseMatrixProductBench<float, tags::Cell> DMPBCF("CELL DenseMatrix Product Benchmark - float", 54);
DenseMatrixProductBench<double, tags::Cell> DMPBCD("CELL DenseMatrix Product Benchmark - double", 20);
#endif


template <typename DT_, typename Tag_>
class BMDVProductBench :
    public Benchmark
{
    private:
        int _count;

    public:
        BMDVProductBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=18)
            {
                cores.push_back(Tag_::name);
                DenseVector<DT_> dv0((j + 1) * 8192, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 8192, DT_(rand()));
                BandedMatrix<DT_> bm((j + 1) * 8192, dv1);
                bm.insert_band(3, dv1.copy());
                bm.insert_band(-3, dv1.copy());
                bm.insert_band(15, dv1.copy());

                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            Product<Tag_>::value(bm, dv0);
                            }
                            );
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//BMDVProductBench<float, tags::CPU> BMDVPBTPF("SingleCore BandedMatrix DenseVector Relax Product Benchmark - float", 20);
//BMDVProductBench<double, tags::CPU> BMDVPBTPD("SingleCore BandedMatrix DenseVector Relax Product Benchmark - double", 20);
//BMDVProductBench<float, tags::CPU::MultiCore> BMDVPBMCF("MultiCore BandedMatrix DenseVector Relax Product Benchmark - float", 20);
//BMDVProductBench<double, tags::CPU::MultiCore> BMDVPBMCD("MultiCore BandedMatrix DenseVector Relax Product Benchmark - double", 20);
#ifdef HONEI_SSE
BMDVProductBench<float, tags::CPU::SSE> BMDVPBSSEF("SSE BandedMatrix DenseVector Relax Product Benchmark - float", 400);
BMDVProductBench<double, tags::CPU::SSE> BMDVPBSSED("SSE BandedMatrix DenseVector Relax Product Benchmark - double", 400);
BMDVProductBench<float, tags::CPU::MultiCore::SSE> BMDVPBMCSSEF("MultiCore SSE BandedMatrix DenseVector Relax Product Benchmark - float", 400);
BMDVProductBench<double, tags::CPU::MultiCore::SSE> BMDVPBMCSSED("MultiCore SSE BandedMatrix DenseVector Relax Product Benchmark - double", 400);
#endif
#elif HONEI_CELL
BMDVProductBench<float, tags::Cell> BMDVPBCF("CELL BandedMatrix DenseVector Relax Product Benchmark - float", 400);
BMDVProductBench<double, tags::Cell> BMDVPBCD("CELL BandedMatrix DenseVector Relax Product Benchmark - double", 200);
#endif


template <typename DT_, typename Tag_>
class BMDVQ1ProductBench :
    public Benchmark
{
    private:
        int _count;

    public:
        BMDVQ1ProductBench(const std::string & id, int count) :
            Benchmark(id)
    {
        register_tag(Tag_::name);
        _plots = true;
        _count = count;
    }

        virtual void run()
        {
            BenchmarkInfo info;
            std::list<BenchmarkInfo> infolist;
            std::list<std::string> cores;
            for (unsigned long j(0) ; j < _count ; j+=3)
            {
                cores.push_back(Tag_::name);
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

                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(
                            for (unsigned long k(0) ; k < 10 ; ++k)
                            {
                            Product<Tag_>::value(bm, dv0);
                            }
                            );
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info * 10);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
//BMDVQ1ProductBench<float, tags::CPU> BMDVQ1PBTPF("SingleCore BandedMatrix DenseVector Q1 Product Benchmark - float", 20);
//BMDVQ1ProductBench<double, tags::CPU> BMDVQ1PBTPD("SingleCore BandedMatrix DenseVector Q1 Product Benchmark - double", 20);
//BMDVQ1ProductBench<float, tags::CPU::MultiCore> BMDVQ1PBMCF("MultiCore BandedMatrix DenseVector Q1 Product Benchmark - float", 20);
//BMDVQ1ProductBench<double, tags::CPU::MultiCore> BMDVQ1PBMCD("MultiCore BandedMatrix DenseVector Q1 Product Benchmark - double", 20);
#ifdef HONEI_SSE
BMDVQ1ProductBench<float, tags::CPU::SSE> BMDVQ1PBSSEF("SSE BandedMatrix DenseVector Q1 Product Benchmark - float", 150);
BMDVQ1ProductBench<double, tags::CPU::SSE> BMDVQ1PBSSED("SSE BandedMatrix DenseVector Q1 Product Benchmark - double", 150);
BMDVQ1ProductBench<float, tags::CPU::MultiCore::SSE> BMDVQ1PBMCSSEF("MultiCore SSE BandedMatrix DenseVector Q1 Product Benchmark - float", 150);
BMDVQ1ProductBench<double, tags::CPU::MultiCore::SSE> BMDVQ1PBMCSSED("MultiCore SSE BandedMatrix DenseVector Q1 Product Benchmark - double", 150);
#endif
#elif HONEI_CELL
BMDVQ1ProductBench<float, tags::Cell> BMDVQ1PBCF("CELL BandedMatrix DenseVector Q1 Product Benchmark - float", 150);
BMDVQ1ProductBench<double, tags::Cell> BMDVQ1PBCD("CELL BandedMatrix DenseVector Q1 Product Benchmark - double", 150);
#endif

