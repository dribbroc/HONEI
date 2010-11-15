/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef ALLBENCH
#include <benchmark/benchmark.cc>

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
        std::string _m_f, _v_f, _i_f;
    public:
        SMELLCGBench(const std::string & id, unsigned long size, unsigned long count, std::string m_file, std::string v_file, std::string i_file) :
            Benchmark(id)
        {
            register_tag(Tag_::name);
            _size = size;
            _count = count;
            _m_f = m_file;
            _v_f = v_file;
            _i_f = i_file;
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
            SparseMatrix<DataType_> tsmatrix2(rows, columns, r, c, data);
            SparseMatrixELL<DataType_> smatrix2(tsmatrix2);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DataType_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DataType_(0)));

            DenseVector<DataType_> diag_inverted(rows, DataType_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DataType_(1)/data[i];
            }

            unsigned long used_iters(0);
            for(unsigned long i(0) ; i < _count ; ++i)
            {
                std::string filename_3(HONEI_SOURCEDIR);
                filename_3 += "/honei/math/testdata/";
                filename_3 += _i_f;
                DenseVector<DataType_> result(VectorIO<io_formats::EXP>::read_vector(filename_3, DataType_(0)));

                if (Tag_::tag_value == tags::tv_gpu_cuda)
                {
                    result.lock(lm_read_and_write, tags::GPU::CUDA::memory_value);
                    result.unlock(lm_read_and_write);
                }
                //DenseVector<DataType_> result_c(result.copy());
                //Defect<Tag_>::value(result, rhs, smatrix2, result_c);

                BENCHMARK(
                        (ConjugateGradients<Tag_, JAC>::value(smatrix2, rhs, result, diag_inverted, 10000ul, used_iters, 1e-8));
#ifdef HONEI_CUDA
                        if (Tag_::tag_value == tags::tv_gpu_cuda)
                            cuda::GPUPool::instance()->flush();
#endif
                        );
            }
            BenchmarkInfo info;
            BenchmarkInfo info_pre;
            info_pre.flops = 2 * non_zeros +  3 * rows;
            info_pre.load = rows * sizeof(DataType_);
            info_pre.store = rows * sizeof(DataType_);
            info.flops = non_zeros * 2 + 13 * rows + 2;
            info.load = (non_zeros*2 + 13 * rows + 3)* sizeof(DataType_);
            info.store = (5 * rows + 3)* sizeof(DataType_);
            evaluate(info * used_iters + info_pre);
        }
};
#ifdef HONEI_SSE
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_0double7("SM ELL 0 L2 CG Benchmark SSE double: " , 0, 1, "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
/*SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_0double7("SM ELL 0 L7 CG Benchmark SSE double: " , 0, 1, "l7/area51_full_0.m", "l7/area51_rhs_0", "l7/area51_init_0");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_0double8("SM ELL 0 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_init_0");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_0double9("SM ELL 0 L9 CG Benchmark SSE double: " , 0, 1, "l9/area51_full_0.m", "l9/area51_rhs_0", "l9/area51_init_0");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_0double10("SM ELL 0 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_0.m", "l10/area51_rhs_0", "l10/area51_init_0");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_1double7("SM ELL 1 L7 CG Benchmark SSE double: " , 0, 1, "l7/area51_full_1.m", "l7/area51_rhs_1", "l7/area51_init_1");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_1double8("SM ELL 1 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_1.m", "l8/area51_rhs_1", "l8/area51_init_1");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_1double9("SM ELL 1 L9 CG Benchmark SSE double: " , 0, 1, "l9/area51_full_1.m", "l9/area51_rhs_1", "l9/area51_init_1");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_1double10("SM ELL 1 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_1.m", "l10/area51_rhs_1", "l10/area51_init_1");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_2double7("SM ELL 2 L7 CG Benchmark SSE double: " , 0, 1, "l7/area51_full_2.m", "l7/area51_rhs_2", "l7/area51_init_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_2double8("SM ELL 2 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2", "l8/area51_init_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_2double9("SM ELL 2 L9 CG Benchmark SSE double: " , 0, 1, "l9/area51_full_2.m", "l9/area51_rhs_2", "l9/area51_init_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_2double10("SM ELL 2 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2", "l10/area51_init_2");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_3double7("SM ELL 3 L7 CG Benchmark SSE double: " , 0, 1, "l7/area51_full_3.m", "l7/area51_rhs_3", "l7/area51_init_3");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_3double8("SM ELL 3 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_3.m", "l8/area51_rhs_3", "l8/area51_init_3");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_3double9("SM ELL 3 L9 CG Benchmark SSE double: " , 0, 1, "l9/area51_full_3.m", "l9/area51_rhs_3", "l9/area51_init_3");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_3double10("SM ELL 3 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_3.m", "l10/area51_rhs_3", "l10/area51_init_3");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_4double7("SM ELL 4 L7 CG Benchmark SSE double: " , 0, 1, "l7/area51_full_4.m", "l7/area51_rhs_4", "l7/area51_init_4");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_4double8("SM ELL 4 L8 CG Benchmark SSE double: " , 0, 1, "l8/area51_full_4.m", "l8/area51_rhs_4", "l8/area51_init_4");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_4double9("SM ELL 4 L9 CG Benchmark SSE double: " , 0, 1, "l9/area51_full_4.m", "l9/area51_rhs_4", "l9/area51_init_4");
SMELLCGBench<tags::CPU::SSE, double> sseSMELLCGBench_4double10("SM ELL 4 L10 CG Benchmark SSE double: " , 0, 1, "l10/area51_full_4.m", "l10/area51_rhs_4", "l10/area51_init_4");*/
/*SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_0double7("SM ELL 0 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_0double7("SM ELL 0 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l7/area51_full_0.m", "l7/area51_rhs_0", "l7/area51_init_0");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_0double8("SM ELL 0 L8 CG Benchmark MultiCore::SSE double: " , 0, 1, "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_init_0");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_0double9("SM ELL 0 L9 CG Benchmark MultiCore::SSE double: " , 0, 1, "l9/area51_full_0.m", "l9/area51_rhs_0", "l9/area51_init_0");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_0double10("SM ELL 0 L10 CG Benchmark MultiCore::SSE double: " , 0, 1, "l10/area51_full_0.m", "l10/area51_rhs_0", "l10/area51_init_0");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_1double7("SM ELL 1 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l7/area51_full_1.m", "l7/area51_rhs_1", "l7/area51_init_1");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_1double8("SM ELL 1 L8 CG Benchmark MultiCore::SSE double: " , 0, 1, "l8/area51_full_1.m", "l8/area51_rhs_1", "l8/area51_init_1");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_1double9("SM ELL 1 L9 CG Benchmark MultiCore::SSE double: " , 0, 1, "l9/area51_full_1.m", "l9/area51_rhs_1", "l9/area51_init_1");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_1double10("SM ELL 1 L10 CG Benchmark MultiCore::SSE double: " , 0, 1, "l10/area51_full_1.m", "l10/area51_rhs_1", "l10/area51_init_1");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_2double7("SM ELL 2 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l7/area51_full_2.m", "l7/area51_rhs_2", "l7/area51_init_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_2double8("SM ELL 2 L8 CG Benchmark MultiCore::SSE double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2", "l8/area51_init_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_2double9("SM ELL 2 L9 CG Benchmark MultiCore::SSE double: " , 0, 1, "l9/area51_full_2.m", "l9/area51_rhs_2", "l9/area51_init_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_2double10("SM ELL 2 L10 CG Benchmark MultiCore::SSE double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2", "l10/area51_init_2");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_3double7("SM ELL 3 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l7/area51_full_3.m", "l7/area51_rhs_3", "l7/area51_init_3");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_3double8("SM ELL 3 L8 CG Benchmark MultiCore::SSE double: " , 0, 1, "l8/area51_full_3.m", "l8/area51_rhs_3", "l8/area51_init_3");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_3double9("SM ELL 3 L9 CG Benchmark MultiCore::SSE double: " , 0, 1, "l9/area51_full_3.m", "l9/area51_rhs_3", "l9/area51_init_3");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_3double10("SM ELL 3 L10 CG Benchmark MultiCore::SSE double: " , 0, 1, "l10/area51_full_3.m", "l10/area51_rhs_3", "l10/area51_init_3");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_4double7("SM ELL 4 L7 CG Benchmark MultiCore::SSE double: " , 0, 1, "l7/area51_full_4.m", "l7/area51_rhs_4", "l7/area51_init_4");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_4double8("SM ELL 4 L8 CG Benchmark MultiCore::SSE double: " , 0, 1, "l8/area51_full_4.m", "l8/area51_rhs_4", "l8/area51_init_4");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_4double9("SM ELL 4 L9 CG Benchmark MultiCore::SSE double: " , 0, 1, "l9/area51_full_4.m", "l9/area51_rhs_4", "l9/area51_init_4");
SMELLCGBench<tags::CPU::MultiCore::SSE, double> msseSMELLCGBench_4double10("SM ELL 4 L10 CG Benchmark MultiCore::SSE double: " , 0, 1, "l10/area51_full_4.m", "l10/area51_rhs_4", "l10/area51_init_4");*/
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_0double7("SM ELL 0 L2 CG Benchmark CUDA double: " , 0, 1, "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
/*SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_0double7("SM ELL 0 L7 CG Benchmark CUDA double: " , 0, 1, "l7/area51_full_0.m", "l7/area51_rhs_0", "l7/area51_init_0");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_0double8("SM ELL 0 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_init_0");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_0double9("SM ELL 0 L9 CG Benchmark CUDA double: " , 0, 1, "l9/area51_full_0.m", "l9/area51_rhs_0", "l9/area51_init_0");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_0double10("SM ELL 0 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_0.m", "l10/area51_rhs_0", "l10/area51_init_0");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_1double7("SM ELL 1 L7 CG Benchmark CUDA double: " , 0, 1, "l7/area51_full_1.m", "l7/area51_rhs_1", "l7/area51_init_1");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_1double8("SM ELL 1 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_1.m", "l8/area51_rhs_1", "l8/area51_init_1");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_1double9("SM ELL 1 L9 CG Benchmark CUDA double: " , 0, 1, "l9/area51_full_1.m", "l9/area51_rhs_1", "l9/area51_init_1");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_1double10("SM ELL 1 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_1.m", "l10/area51_rhs_1", "l10/area51_init_1");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_2double7("SM ELL 2 L7 CG Benchmark CUDA double: " , 0, 1, "l7/area51_full_2.m", "l7/area51_rhs_2", "l7/area51_init_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_2double8("SM ELL 2 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_2.m", "l8/area51_rhs_2", "l8/area51_init_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_2double9("SM ELL 2 L9 CG Benchmark CUDA double: " , 0, 1, "l9/area51_full_2.m", "l9/area51_rhs_2", "l9/area51_init_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_2double10("SM ELL 2 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_2.m", "l10/area51_rhs_2", "l10/area51_init_2");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_3double7("SM ELL 3 L7 CG Benchmark CUDA double: " , 0, 1, "l7/area51_full_3.m", "l7/area51_rhs_3", "l7/area51_init_3");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_3double8("SM ELL 3 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_3.m", "l8/area51_rhs_3", "l8/area51_init_3");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_3double9("SM ELL 3 L9 CG Benchmark CUDA double: " , 0, 1, "l9/area51_full_3.m", "l9/area51_rhs_3", "l9/area51_init_3");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_3double10("SM ELL 3 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_3.m", "l10/area51_rhs_3", "l10/area51_init_3");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_4double7("SM ELL 4 L7 CG Benchmark CUDA double: " , 0, 1, "l7/area51_full_4.m", "l7/area51_rhs_4", "l7/area51_init_4");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_4double8("SM ELL 4 L8 CG Benchmark CUDA double: " , 0, 1, "l8/area51_full_4.m", "l8/area51_rhs_4", "l8/area51_init_4");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_4double9("SM ELL 4 L9 CG Benchmark CUDA double: " , 0, 1, "l9/area51_full_4.m", "l9/area51_rhs_4", "l9/area51_init_4");
SMELLCGBench<tags::GPU::CUDA, double> cudaSMELLCGBench_4double10("SM ELL 4 L10 CG Benchmark CUDA double: " , 0, 1, "l10/area51_full_4.m", "l10/area51_rhs_4", "l10/area51_init_4");*/
#endif
#endif
