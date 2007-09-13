/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <string>
#include <tr1/memory>
#endif

#include <libla/dot_product.hh>

using namespace std;
using namespace honei;


template <typename DataType_>

class DotProductBench :
    public Benchmark
{
    private:
        unsigned long _size;
        int _count;
    public:
        DotProductBench(const std::string & id, unsigned long size, int count) :
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
                DenseVector<DataType_> dv0(_size, DataType_(rand()));
                DenseVector<DataType_> dv1(_size, DataType_(rand()));
                BENCHMARK(p0 = DotProduct<DataType_>::value(dv1,dv0));
            }
            cout << endl;
            evaluate(2*_size);
        }
};

// ScalarProductBench<float> SPBenchfloat1("Scalar Product Benchmark - vector size: 1,000,000, float", 1000000, 10);
// ScalarProductBench<double> SPBenchdouble1("Scalar Product Benchmark - vector size: 1,000,000, double", 1000000, 10);
DotProductBench<float> DPBenchfloat("Dot Product Benchmark - vector size: 10,000 float", 10000, 10);
DotProductBench<double> DPBenchdouble("Dot Product Benchmark - vector size: 10,000 double", 10000, 10);
