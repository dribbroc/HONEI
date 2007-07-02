/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.cc>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_scaled_sum.hh>
#include <libla/vector_norm.hh>

#include <tr1/memory>
#include <string>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class DenseVectorScaledSumBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseVectorScaledSumBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv0(new DenseVector<DataType_>(_size, static_cast<DataType_>(rand())));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(_size, static_cast<DataType_>(rand())));
                DataType_ alpha(static_cast<DataType_>(rand()));
                DataType_ beta(static_cast<DataType_>(rand()));
                BENCHMARK(DenseVector<DataType_> sum1(VectorScaledSum<>::value(*dv0, *dv1, alpha, beta)););
            }
            evaluate(3*_size);
        }
};

DenseVectorScaledSumBench<float> DVSSBenchfloat1("Dense Vector Scaled Sum Benchmark - vector size: 1,000,000, float", 1000000, 10);
DenseVectorScaledSumBench<double> DVSSBenchdouble1("Dense Vector Scaled Sum Benchmark - vector size: 1,000,000, double", 1000000, 10);
DenseVectorScaledSumBench<float> DVSSBenchfloat2("Dense Vector Scaled Sum Benchmark - vector size: 10,000,000, float", 10000000, 10);
DenseVectorScaledSumBench<double> DVSSBenchdouble2("Dense Vector Scaled Sum Benchmark - vector size: 10,000,000, double", 10000000, 10);
