/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.cc>
#include <tr1/memory>
#include <string>
#endif

#include <libla/matrix_product.hh>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class DenseMatrixProductBench :
    public Benchmark
{
    private:
        int _size;
        int _count;
    public:
        DenseMatrixProductBench(const std::string & id, int size, int count) :
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
                BENCHMARK(MatrixProduct<DataType_>::value(*dm0, *dm1));
            }
            evaluate(_size*_size*_size*2);
        }
};

DenseMatrixProductBench<float> DMPBenchfloat1("Dense Matrix Product Benchmark - matrix size: 128x128, float", 128, 10);
DenseMatrixProductBench<double> DMPBenchdouble1("Dense Matrix Product Benchmark - matrix size: 128x128, double", 128, 10);
DenseMatrixProductBench<float> DMPBenchfloat2("Dense Matrix Product Benchmark - matrix size: 256x256, float", 256, 10);
DenseMatrixProductBench<double> DMPBenchdouble2("Dense Matrix Product Benchmark - matrix size: 256x256, double", 256, 10);
