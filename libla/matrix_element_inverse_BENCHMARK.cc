/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.cc>
#include <libla/dense_matrix.hh>
#include <libla/matrix_element_inverse.hh>
#include <libla/matrix_error.cc>

#include <tr1/memory>
#include <string>
 
using namespace std;
using namespace pg512;
 

template <typename DataType_>

class DenseMatrixElementInverseBench :
    public Benchmark
{
    private:
        int _size;
        int _count;

    public:
        DenseMatrixElementInverseBench(const std::string & id, int size, int count) :
            Benchmark(id)
        {
            _size  = size;
            _count = count;
        }

        virtual void run()
        {
            for(int i = 0; i < _count; ++i)
        {

            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm0(new DenseMatrix<DataType_>(_size, _size, DataType_(12)) );
            BENCHMARK(MatrixElementInverse<DataType_>::value(*dm0));
        }
        evaluate(_size * _size);
    }
};

DenseMatrixElementInverseBench<float>  MEIBenchfloat ("Matrix Element Inverse Benchmark: size: 1000x1000, float",  1000, 10);
DenseMatrixElementInverseBench<double> MEIBenchdouble("Matrix Element Inverse Benchmark: size: 1000x1000, double", 1000, 10);
