/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>
#include <tr1/memory>
#include <string>
#endif

#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/jacobi.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/operations.hh>
#include <iostream>
#include <cmath>
//using namespace std;
using namespace honei;

template <typename Tag_, typename DataType_>
class SMELLJacobiBench :
    public Benchmark
{
    private:
        unsigned long _size;
        unsigned long _count;
        std::string _m_f, _v_f;
    public:
        SMELLJacobiBench(const std::string & id, unsigned long size, unsigned long count, std::string m_file, std::string v_file) :
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

            //std::string filename(HONEI_SOURCEDIR);
            //filename += "/honei/math/testdata/5pt_10x10.mtx";
            //filename += "/honei/math/testdata/test_0.mtx";
            std::string filename = "/home/user/mgeveler/nobackup/feat2/Featflow2/area51/renumbenchmark/";
            filename += _m_f;
            unsigned long non_zeros(MatrixIO<io_formats::M>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DataType_> data(non_zeros);

            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrixELL<DataType_> smatrix2(rows, columns, r, c, data);

            std::string filename_2 = "/home/user/mgeveler/nobackup/feat2/Featflow2/area51/renumbenchmark/";
            filename_2 += _v_f;
            DenseVector<DataType_> rhs(rows, DataType_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            DenseVector<DataType_> initial_guess(rows, DataType_(1));
            DenseVector<DataType_> diag_inverted(rows, DataType_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DataType_(1)/data[i];
            }

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                BENCHMARK(
                        Jacobi<Tag_>::value(initial_guess, smatrix2, rhs, 100ul, DataType_(0.7), diag_inverted);
                        );
#ifdef HONEI_CUDA
                        cuda_thread_synchronize();
#endif
            }
            BenchmarkInfo info;
            info.flops = non_zeros * 2 + 3 * rows;
            evaluate(info * 100);
        }
};
/*
#ifdef HONEI_CUDA
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda0("SM ELL 0 Jacobi Benchmark CUDA float: " , 0, 10, "l10/area51_full_0.m", "l10/area51_rhs_0" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda1("SM ELL 1 Jacobi Benchmark CUDA float: " , 0, 10, "l10/area51_full_1.m", "l10/area51_rhs_1" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda2("SM ELL 2 Jacobi Benchmark CUDA float: " , 0, 10, "l10/area51_full_2.m", "l10/area51_rhs_2" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda3("SM ELL 3 Jacobi Benchmark CUDA float: " , 0, 10, "l10/area51_full_3.m", "l10/area51_rhs_3" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda4("SM ELL 4 Jacobi Benchmark CUDA float: " , 0, 10, "l10/area51_full_4.m", "l10/area51_rhs_4" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda5("SM ELL 0f Jacobi Benchmark CUDA float: " , 0, 10, "l10_full/area51_full_0.m", "l10_full/area51_rhs_0" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda6("SM ELL 1f Jacobi Benchmark CUDA float: " , 0, 10, "l10_full/area51_full_1.m", "l10_full/area51_rhs_1" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda7("SM ELL 2f Jacobi Benchmark CUDA float: " , 0, 10, "l10_full/area51_full_2.m", "l10_full/area51_rhs_2" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda8("SM ELL 3f Jacobi Benchmark CUDA float: " , 0, 10, "l10_full/area51_full_3.m", "l10_full/area51_rhs_3" );
SMELLJacobiBench<tags::GPU::CUDA, float> SMELLDVPBench_float_cuda9("SM ELL 4f Jacobi Benchmark CUDA float: " , 0, 10, "l10_full/area51_full_4.m", "l10_full/area51_rhs_4" );
#ifdef HONEI_CUDA_DOUBLE
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double0("SM ELL 0 Jacobi Benchmark CUDA double: " , 0, 10, "l10/area51_full_0.m", "l10/area51_rhs_0" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double1("SM ELL 1 Jacobi Benchmark CUDA double: " , 0, 10, "l10/area51_full_1.m", "l10/area51_rhs_1" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double2("SM ELL 2 Jacobi Benchmark CUDA double: " , 0, 10, "l10/area51_full_2.m", "l10/area51_rhs_2" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double3("SM ELL 3 Jacobi Benchmark CUDA double: " , 0, 10, "l10/area51_full_3.m", "l10/area51_rhs_3" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double4("SM ELL 4 Jacobi Benchmark CUDA double: " , 0, 10, "l10/area51_full_4.m", "l10/area51_rhs_4" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double5("SM ELL 0f Jacobi Benchmark CUDA double: " , 0, 10, "l10_full/area51_full_0.m", "l10_full/area51_rhs_0" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double6("SM ELL 1f Jacobi Benchmark CUDA double: " , 0, 10, "l10_full/area51_full_1.m", "l10_full/area51_rhs_1" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double7("SM ELL 2f Jacobi Benchmark CUDA double: " , 0, 10, "l10_full/area51_full_2.m", "l10_full/area51_rhs_2" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double8("SM ELL 3f Jacobi Benchmark CUDA double: " , 0, 10, "l10_full/area51_full_3.m", "l10_full/area51_rhs_3" );
SMELLJacobiBench<tags::GPU::CUDA, double> SMELLDVPBench_double9("SM ELL 4f Jacobi Benchmark CUDA double: " , 0, 10, "l10_full/area51_full_4.m", "l10_full/area51_rhs_4" );
#endif
#endif
*/
