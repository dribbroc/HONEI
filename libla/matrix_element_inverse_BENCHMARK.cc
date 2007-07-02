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

class MatrixElementInverseBench :
    public Benchmark
{
    private:
        int _size;
        int _count;

    public:
        MatrixElementInverseBench(const std::string & id, int size, int count) :
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

MatrixElementInverseBench<float>  MEIBenchfloat ("Matrix Element Inverse Benchmark: size: 2000x2000, float",  5000, 10);
MatrixElementInverseBench<double> MEIBenchdouble("Matrix Element Inverse Benchmark: size: 2000x2000, double", 5000, 10);
