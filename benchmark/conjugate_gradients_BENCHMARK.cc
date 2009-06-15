/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/conjugate_gradients.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class SMELLCGBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
        std::string _m_f, _v_f;
    public:
        SMELLCGBench(const std::string & id, unsigned long size, unsigned long count, std::string m_file, std::string v_file) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
            _m_f = m_file;
            _v_f = v_file;
        }

        virtual void run()
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            //std::string filename = "/home/mgeveler/testdata/";
            filename += _m_f;
            unsigned long non_zeros(MatrixIO<io_formats::M>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DataType_> data(non_zeros);

            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrixELL<DataType_> smatrix2(rows, columns, r, c, data);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DataType_> rhs(rows, DataType_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            DenseVector<DataType_> diag_inverted(rows, DataType_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DataType_(1)/data[i];
            }

            DenseVector<DataType_> result(rhs.size(), DataType_(1));

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        (ConjugateGradients<Tag_, JAC>::value(smatrix2, rhs, result, diag_inverted, 100ul));
                        );
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
            }
            BenchmarkInfo info;
            BenchmarkInfo info_pre;
            info_pre.flops = ((2 * non_zeros + 6 * rows));
            info.flops = non_zeros * 2 + 13 * rows + 2;
            evaluate(info * 100 + info_pre);
        }
};
#ifdef HONEI_CUDA
SMELLCGBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda2("SM ELL 2 CG Benchmark CUDA float: " , 0, 1, "area51_full_0.m", "area51_rhs_0");
#ifdef HONEI_CUDA_DOUBLE
SMELLCGBench<tags::GPU::CUDA, double> SMELLDVPBench_double_cuda2("SM ELL 2 CG Benchmark CUDA double: " , 0, 1, "area51_full_0.m", "area51_rhs_0");
#endif
#endif
