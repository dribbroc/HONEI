/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <libla/dense_vector.hh>
#include <tr1/memory>
#include <string>
#endif

#include <libla/scalar_product.hh>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class ScalarProductBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        ScalarProductBench(const std::string & id, int size, int count) :
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
                std::tr1::shared_ptr<DenseVector<DataType_> > dv0(new DenseVector<DataType_>(_size, static_cast<DataType_>(rand())));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(_size, static_cast<DataType_>(rand())));
                BENCHMARK(p0 = ScalarProduct<DataType_>::value(*dv1,*dv0));
            }
            evaluate(2*_size);
        }
};

ScalarProductBench<float> SPBenchfloat1("Scalar Product Benchmark - vector size: 1,000,000, float", 1000000, 10);
ScalarProductBench<double> SPBenchdouble1("Scalar Product Benchmark - vector size: 1,000,000, double", 1000000, 10);
ScalarProductBench<float> SPBenchfloat2("Scalar Product Benchmark - vector size: 10,000,000, float", 10000000, 10);
ScalarProductBench<double> SPBenchdouble2("Scalar Product Benchmark - vector size: 10,000,000, double", 10000000, 10);
