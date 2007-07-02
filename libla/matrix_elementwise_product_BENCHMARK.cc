/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_elementwise_product.hh>
#include <libla/matrix_error.cc>

#include <tr1/memory>
#include <string>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class DenseMatrixElementwiseProductBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseMatrixElementwiseProductBench(const std::string & id, int size, int count) :
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
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm0(new DenseMatrix<DataType_>(_size, _size, DataType_(rand())) );
                std::tr1::shared_ptr<DenseMatrix<DataType_> > dm1(new DenseMatrix<DataType_>(_size, _size, DataType_(rand())) );
                BENCHMARK(MatrixElementwiseProduct<DataType_>::value(*dm0, *dm1));
            }
            evaluate(_size*_size);
        }
};

DenseMatrixElementwiseProductBench<float> DMEPBenchfloat1("Dense Matrix Elementwise Product Benchmark - matrix size: 2048x2048, float", 2048, 10);
DenseMatrixElementwiseProductBench<double> DMEPBenchdouble1("Dense Matrix Elementwise Product Benchmark - matrix size: 2048x2048, double", 2048, 10);
DenseMatrixElementwiseProductBench<float> DMEPBenchfloat2("Dense Matrix Elementwise Product Benchmark - matrix size: 4096x4096, float", 4096, 10);
DenseMatrixElementwiseProductBench<double> DMEPBenchdouble2("Dense Matrix Elementwise Product Benchmark - matrix size: 4096x4096, double", 4096, 10);
