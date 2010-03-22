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
#include <honei/backends/cuda/gpu_pool.hh>
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

            for(unsigned long i(0) ; i < _count ; ++i)
            {
                DenseVector<DataType_> result(rhs.size(), DataType_(1));
                if (Tag_::tag_value == tags::tv_gpu_cuda)
                {
                    result.lock(lm_read_and_write, tags::GPU::CUDA::memory_value);
                    result.unlock(lm_read_and_write);
                }
                BENCHMARK(
                        (ConjugateGradients<Tag_, JAC>::value(smatrix2, rhs, result, diag_inverted, 1000ul));
                        );
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
            }
            BenchmarkInfo info;
            BenchmarkInfo info_pre;
            info_pre.flops = 2 * non_zeros +  3 * rows;
            info_pre.load = rows * sizeof(DataType_);
            info_pre.store = rows * sizeof(DataType_);
            info.flops = non_zeros * 2 + 13 * rows + 2;
            info.load = (non_zeros*2 + 13 * rows + 3)* sizeof(DataType_);
            info.store = (5 * rows + 3)* sizeof(DataType_);
            evaluate(info * 1000 + info_pre);
        }
};
#ifdef HONEI_SSE
SMELLCGBench<tags::CPU::SSE, float> sseSMELLCGBench_float2("SM ELL 2 L2 CG Benchmark SSE float: " , 0, 1, "l2/area51_full_2.m", "l2/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, float> sseSMELLCGBench_float8("SM ELL 2 L8 CG Benchmark SSE float: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, float> sseSMELLCGBench_float10("SM ELL 2 L10 CG Benchmark SSE float: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, float> sseSMELLCGBench_float11("SM ELL 2 L11 CG Benchmark SSE float: " , 0, 1, "l11/area51_full_2.m", "l11/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_double2("SM ELL 2 L2 CG Benchmark SSE double: " , 0, 1, "l2/area51_full_2.m", "l2/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_double8("SM ELL 2 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_double10("SM ELL 2 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_double11("SM ELL 2 L11 CG Benchmark SSE double: " , 0, 1, "l11/area51_full_2.m", "l11/area51_rhs_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> mcsseSMELLCGBench_double8("MC SM ELL 2 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> mcsseSMELLCGBench_double10("MC SM ELL 2 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2");
#endif
#ifdef HONEI_CUDA
SMELLCGBench<tags::GPU::CUDA, float> cudaSMELLCGBench_float2("SM ELL 2 L2 CG Benchmark CUDA float: " , 0, 1, "l2/area51_full_2.m", "l2/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, float> cudaSMELLCGBench_float8("SM ELL 2 L8 CG Benchmark CUDA float: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, float> cudaSMELLCGBench_float10("SM ELL 2 L10 CG Benchmark CUDA float: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, float> cudaSMELLCGBench_float11("SM ELL 2 L11 CG Benchmark CUDA float: " , 0, 1, "l11/area51_full_2.m", "l11/area51_rhs_2");
#ifdef HONEI_CUDA_DOUBLE
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_double2("SM ELL 2 L2 CG Benchmark CUDA double: " , 0, 1, "l2/area51_full_2.m", "l2/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_double8("SM ELL 2 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_double10("SM ELL 2 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_double11("SM ELL 2 L11 CG Benchmark CUDA double: " , 0, 1, "l11/area51_full_2.m", "l11/area51_rhs_2");
#endif
#endif
