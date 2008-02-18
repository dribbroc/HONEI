#include <benchmark/benchmark.cc>

#include <honei/libla/sum.hh>
#include <honei/libla/product.hh>
#include <honei/libla/reduction.hh>
#include <honei/libla/dot_product.hh>



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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv((j + 1) * 131072, DT_(rand()));
                DT_ dt;
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK((dt = Reduction<rt_sum, Tag_>::value(dv)));
                }
                info = Reduction<rt_sum>::get_benchmark_info(dv);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
DenseVectorRTSBench<float, tags::CPU> DVRBF("SingleCore DenseVector Reduction to Sum Benchmark - float", 20);
DenseVectorRTSBench<double, tags::CPU> DVRBD("SingleCore DenseVector Reduction to Sum Benchmark - double", 20);
DenseVectorRTSBench<float, tags::CPU::MultiCore> DVRBMCF("MultiCore DenseVector Reduction to Sum Benchmark - float", 20);
DenseVectorRTSBench<double, tags::CPU::MultiCore> DVRBMCD("MultiCore DenseVector Reduction to Sum Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorRTSBench<float, tags::CPU::SSE> DVRBSSEF("SSE DenseVector Reduction to Sum Benchmark - float", 40);
DenseVectorRTSBench<double, tags::CPU::SSE> DVRBSSED("SSE DenseVector Reduction to Sum Benchmark - double", 40);
DenseVectorRTSBench<float, tags::CPU::MultiCore::SSE> DVRBMCSSEF("MultiCore SSE DenseVector Reduction to Sum Benchmark - float", 40);
DenseVectorRTSBench<double, tags::CPU::MultiCore::SSE> DVRBMCSSED("MultiCore SSE DenseVector Reduction to Sum Benchmark - double", 40);
#endif
#elif HONEI_CELL
DenseVectorRTSBench<float, tags::Cell> DVRBCF("CELL DenseVector Reduction to Sum Benchmark - float", 40);
DenseVectorRTSBench<double, tags::Cell> DVRBCD("CELL DenseVector Reduction to Sum Benchmark - double", 40);
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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(DotProduct<Tag_>::value(dv0, dv1));
                }
                info = DotProduct<>::get_benchmark_info(dv0, dv1);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
DenseVectorDotProductBench<float, tags::CPU> DVDPBF("SingleCore DenseVector DotProduct Benchmark - float", 20);
DenseVectorDotProductBench<double, tags::CPU> DVDPBD("SingleCore DenseVector DotProduct Benchmark - double", 20);
DenseVectorDotProductBench<float, tags::CPU::MultiCore> DVDPBMCF("MultiCore DenseVector DotProduct Benchmark - float", 20);
DenseVectorDotProductBench<double, tags::CPU::MultiCore> DVDPBMCD("MultiCore DenseVector DotProduct Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorDotProductBench<float, tags::CPU::SSE> DVDPBSSEF("SSE DenseVector DotProduct Benchmark - float", 40);
DenseVectorDotProductBench<double, tags::CPU::SSE> DVDPBSSED("SSE DenseVector DotProduct Benchmark - double", 40);
DenseVectorDotProductBench<float, tags::CPU::MultiCore::SSE> DVDPBMCSSEF("MultiCore SSE DenseVector DotProduct Benchmark - float", 40);
DenseVectorDotProductBench<double, tags::CPU::MultiCore::SSE> DVDPBMCSSED("MultiCore SSE DenseVector DotProduct Benchmark - double", 40);
#endif
#elif HONEI_CELL
DenseVectorDotProductBench<float, tags::Cell> DVDPBCF("CELL DenseVector DotProduct Benchmark - float", 40);
DenseVectorDotProductBench<double, tags::Cell> DVDPBCD("CELL DenseVector DotProduct Benchmark - double", 40);
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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv0((j + 1) * 131072, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 131072, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(Sum<Tag_>::value(dv0, dv1));
                }
                info = Sum<>::get_benchmark_info(dv0, dv1);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
DenseVectorSumBench<float, tags::CPU> DVSBF("SingleCore DenseVector Sum Benchmark - float", 20);
DenseVectorSumBench<double, tags::CPU> DVSBD("SingleCore DenseVector Sum Benchmark - double", 20);
DenseVectorSumBench<float, tags::CPU::MultiCore> DVSBMCF("MultiCore DenseVector Sum Benchmark - float", 20);
DenseVectorSumBench<double, tags::CPU::MultiCore> DVSBMCD("MultiCore DenseVector Sum Benchmark - double", 20);
#ifdef HONEI_SSE
DenseVectorSumBench<float, tags::CPU::SSE> DVSBSSEF("SSE DenseVector Sum Benchmark - float", 40);
DenseVectorSumBench<double, tags::CPU::SSE> DVSBSSED("SSE DenseVector Sum Benchmark - double", 40);
DenseVectorSumBench<float, tags::CPU::MultiCore::SSE> DVSBMCSSEF("MultiCore SSE DenseVector Sum Benchmark - float", 40);
DenseVectorSumBench<double, tags::CPU::MultiCore::SSE> DVSBMCSSED("MultiCore SSE DenseVector Sum Benchmark - double", 40);
#endif
#elif HONEI_CELL
DenseVectorSumBench<float, tags::Cell> DVSBCF("CELL DenseVector Sum Benchmark - float", 40);
DenseVectorSumBench<double, tags::Cell> DVSBCD("CELL DenseVector Sum Benchmark - double", 40);
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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv((j + 1) * 128, DT_(rand()));
                DenseMatrix<DT_> dm((j + 1) * 128, (j + 1) * 128, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(Product<Tag_>::value(dm, dv));
                }
                info = Product<>::get_benchmark_info(dm, dv);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
DMDVProductBench<float, tags::CPU> DMDVPBTPF("SingleCore DenseMatrix DenseVector Product Benchmark - float", 20);
DMDVProductBench<double, tags::CPU> DMDVPBTPD("SingleCore DenseMatrix DenseVector Product Benchmark - double", 20);
DMDVProductBench<float, tags::CPU::MultiCore> DMDVPBMCF("MultiCore DenseMatrix DenseVector Product Benchmark - float", 20);
DMDVProductBench<double, tags::CPU::MultiCore> DMDVPBMCD("MultiCore DenseMatrix DenseVector Product Benchmark - double", 20);
#ifdef HONEI_SSE
DMDVProductBench<float, tags::CPU::SSE> DMDVPBSSEF("SSE DenseMatrix DenseVector Product Benchmark - float", 40);
DMDVProductBench<double, tags::CPU::SSE> DMDVPBSSED("SSE DenseMatrix DenseVector Product Benchmark - double", 40);
DMDVProductBench<float, tags::CPU::MultiCore::SSE> DMDVPBMCSSEF("MultiCore SSE DenseMatrix DenseVector Product Benchmark - float", 40);
DMDVProductBench<double, tags::CPU::MultiCore::SSE> DMDVPBMCSSED("MultiCore SSE DenseMatrix DenseVector Product Benchmark - double", 40);
#endif
#elif HONEI_CELL
DMDVProductBench<float, tags::Cell> DMDVPBCF("CELL DenseMatrix DenseVector Product Benchmark - float", 40);
//DMDVProductBench<double, tags::Cell> DMDVPBCD("CELL DenseMatrix DenseVector Product Benchmark - double", 40);
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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseMatrix<DT_> dm0((j + 1) * 16, (j + 1) * 16, DT_(rand()));
                DenseMatrix<DT_> dm1((j + 1) * 16, (j + 1) * 16, DT_(rand()));
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(Product<Tag_>::value(dm0, dm1));
                }
                info = Product<>::get_benchmark_info(dm0, dm1);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
DenseMatrixProductBench<float, tags::CPU> DMPBTPF("SingleCore DenseMatrix Product Benchmark - float", 15);
DenseMatrixProductBench<double, tags::CPU> DMPBTPD("SingleCore DenseMatrix Product Benchmark - double", 15);
DenseMatrixProductBench<float, tags::CPU::MultiCore> DMPBMCF("MultiCore DenseMatrix Product Benchmark - float", 15);
DenseMatrixProductBench<double, tags::CPU::MultiCore> DMPBMCD("MultiCore DenseMatrix Product Benchmark - double", 15);
#ifdef HONEI_SSE
DenseMatrixProductBench<float, tags::CPU::SSE> DMPBSSEF("SSE DenseMatrix Product Benchmark - float", 40);
DenseMatrixProductBench<double, tags::CPU::SSE> DMPBSSED("SSE DenseMatrix Product Benchmark - double", 40);
DenseMatrixProductBench<float, tags::CPU::MultiCore::SSE> DMPBMCSSEF("MultiCore SSE DenseMatrix Product Benchmark - float", 40);
DenseMatrixProductBench<double, tags::CPU::MultiCore::SSE> DMPBMCSSED("MultiCore SSE DenseMatrix Product Benchmark - double", 40);
#endif
#elif HONEI_CELL
DenseMatrixProductBench<float, tags::Cell> DMPBCF("CELL DenseMatrix Product Benchmark - float", 40);
//DenseMatrixProductBench<double, tags::Cell> DMPBCD("CELL DenseMatrix Product Benchmark - double", 40);
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
            std::list<int> cores;
            for (unsigned long j(0) ; j < _count ; ++j)
            {
                cores.push_back(1);
                DenseVector<DT_> dv0((j + 1) * 8192, DT_(rand()));
                DenseVector<DT_> dv1((j + 1) * 8192, DT_(rand()));
                BandedMatrix<DT_> bm((j + 1) * 8192, dv1);
                for (int i = 1; i < 14 ; i++)
                {
                    bm.insert_band(i * 3, dv1.copy());
                    bm.insert_band(-1 * 5 * i, dv1.copy());
                }
                for(int i(0) ; i < 20 ; ++i)
                {
                    BENCHMARK(Product<Tag_>::value(bm, dv0));
                }
                info = Product<>::get_benchmark_info(bm, dv0);
                infolist.push_back(info);
                std::cout << "finished run " << j + 1 << " / " << _count << std::endl;
            }
            evaluate_to_plotfile(infolist, cores, 20);
        }
};
#ifndef HONEI_CELL
BMDVProductBench<float, tags::CPU> BMDVPBTPF("SingleCore BandedMatrix DenseVector Product Benchmark - float", 20);
BMDVProductBench<double, tags::CPU> BMDVPBTPD("SingleCore BandedMatrix DenseVector Product Benchmark - double", 20);
BMDVProductBench<float, tags::CPU::MultiCore> BMDVPBMCF("MultiCore BandedMatrix DenseVector Product Benchmark - float", 20);
BMDVProductBench<double, tags::CPU::MultiCore> BMDVPBMCD("MultiCore BandedMatrix DenseVector Product Benchmark - double", 20);
#ifdef HONEI_SSE
BMDVProductBench<float, tags::CPU::SSE> BMDVPBSSEF("SSE BandedMatrix DenseVector Product Benchmark - float", 40);
BMDVProductBench<double, tags::CPU::SSE> BMDVPBSSED("SSE BandedMatrix DenseVector Product Benchmark - double", 40);
BMDVProductBench<float, tags::CPU::MultiCore::SSE> BMDVPBMCSSEF("MultiCore SSE BandedMatrix DenseVector Product Benchmark - float", 40);
BMDVProductBench<double, tags::CPU::MultiCore::SSE> BMDVPBMCSSED("MultiCore SSE BandedMatrix DenseVector Product Benchmark - double", 40);
#endif
#elif HONEI_CELL
BMDVProductBench<float, tags::Cell> BMDVPBCF("CELL BandedMatrix DenseVector Product Benchmark - float", 40);
BMDVProductBench<double, tags::Cell> BMDVPBCD("CELL BandedMatrix DenseVector Product Benchmark - double", 40);
#endif
