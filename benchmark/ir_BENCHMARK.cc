/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

//#define SOLVER_VERBOSE 1
#include <honei/math/ir.hh>
#include <honei/math/methods.hh>
#include <honei/math/restriction.hh>
#include <honei/math/prolongation.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <benchmark/benchmark.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/backends/cuda/gpu_pool.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace std;


template <typename OuterTag_, typename InnerTag_, typename DTOuter_, typename DTInner_, typename CycleType_, typename PreconContType_>
class IRBench:
    public Benchmark
{
    private:
        unsigned long _element_type;
        unsigned long _sorting;
        unsigned long _levels;
        unsigned long _max_iters_inner_mg;
        std::string _precon;
        double _damping;

    public:
        IRBench(const std::string & tag, unsigned long et, unsigned long s, unsigned long l, std::string p, double d, unsigned long m) :
            Benchmark(tag)
        {
            register_tag(OuterTag_::name);
            _element_type = et;
            _sorting = s;
            _levels = l;
            _precon = p;
            _damping = d;
            _max_iters_inner_mg = m;
        }

        virtual void run()
        {
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson_advanced/";

            if(_element_type == 2)
                file += "q2_";

            file += "sort_";

            file += stringify(_sorting);

            file += "/";

            if(typeid(OuterTag_) == typeid(tags::CPU::MultiCore::SSE))
            {
                if(_element_type == 1)
                {
                    Configuration::instance()->set_value("ell::threads", 2);
                }
                else
                {
                    Configuration::instance()->set_value("ell::threads", 16);
                }
            }
            else
            {
                if(_element_type == 1)
                {
                    Configuration::instance()->set_value("ell::threads", 1);
                    Configuration::instance()->set_value("cuda::product_smell_dv_double", 128);
                }
                else
                {
                    Configuration::instance()->set_value("ell::threads", 2);
                    Configuration::instance()->set_value("cuda::product_smell_dv_double", 256);
                }
            }

            std::string a_file;
            a_file += file;
            a_file += "A_";
            a_file += stringify(_levels);
            a_file += ".ell";
            SparseMatrixELL<DTOuter_> A(MatrixIO<io_formats::ELL>::read_matrix(a_file, DTOuter_(0)));

            std::string b_file;
            b_file += file;
            b_file += "rhs_";
            b_file += stringify(_levels);
            DenseVector<DTOuter_> b(VectorIO<io_formats::EXP>::read_vector(b_file, DTOuter_(0)));

            std::string x_file;
            x_file += file;
            x_file += "init_";
            x_file += stringify(_levels);
            DenseVector<DTOuter_> x(VectorIO<io_formats::EXP>::read_vector(x_file, DTOuter_(0)));


            if(typeid(InnerTag_) == typeid(tags::CPU::MultiCore::SSE))
            {
                if(_element_type == 1)
                {
                    Configuration::instance()->set_value("ell::threads", 2);
                }
                else
                {
                    Configuration::instance()->set_value("ell::threads", 16);
                }
            }
            else
            {
                if(_element_type == 1)
                {
                    Configuration::instance()->set_value("ell::threads", 1);
                    Configuration::instance()->set_value("cuda::product_smell_dv_double", 128);
                }
                else
                {
                    Configuration::instance()->set_value("ell::threads", 2);
                    Configuration::instance()->set_value("cuda::product_smell_dv_double", 256);
                }
            }
            MGData<SparseMatrixELL<DTInner_>, DenseVector<DTInner_>, SparseMatrixELL<DTInner_>, PreconContType_, DTInner_>  data(MGUtil<InnerTag_,
                    SparseMatrixELL<DTInner_>,
                    DenseVector<DTInner_>,
                    SparseMatrixELL<DTInner_>,
                    PreconContType_,
                    io_formats::ELL,
                    io_formats::EXP,
                    DTInner_,
                    RISmoother<InnerTag_> >::load_data(file, _levels, _damping, _precon));

            MGUtil<InnerTag_,
                SparseMatrixELL<DTInner_>,
                DenseVector<DTInner_>,
                SparseMatrixELL<DTInner_>,
                PreconContType_,
                io_formats::ELL,
                io_formats::EXP,
                DTInner_,
                RISmoother<InnerTag_> >::configure(data, _max_iters_inner_mg, 1, 4, 4, 1, DTInner_(1e-8));

            OperatorList ol(
                    MGCycleCreation<InnerTag_,
                    CycleType_,
                    CG<InnerTag_, methods::VAR>,
                    RISmoother<InnerTag_>,
                    Restriction<InnerTag_, methods::PROLMAT>,
                    Prolongation<InnerTag_, methods::PROLMAT>,
                    DTInner_>::value(data)
                    );

            unsigned long used(0);
            BENCHMARK(
                    (IRSolver<OuterTag_, InnerTag_, Norm<vnt_l_two, true, OuterTag_> >::value(A, b, x, data, ol, double(1e-8), 100ul, used));
#ifdef HONEI_CUDA
                    if (InnerTag_::tag_value == tags::tv_gpu_cuda || OuterTag_::tag_value == tags::tv_gpu_cuda)
                    cuda::GPUPool::instance()->flush();
#endif
                    );

            std::cout << used << std::endl;
            std::cout << data.used_iters << std::endl;
            std::cout << data.used_iters_coarse << std::endl;

            evaluate();

        }
};

#ifdef HONEI_SSE
#ifdef HONEI_CUDA
//JAC
//q1
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort0_l7_jac_v("IRBench cuda/cuda | q1 | sort 0 | L7 | jac | V", 1 , 0, 7, "jac", 0.7, 8);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort0_l8_jac_v("IRBench cuda/cuda | q1 | sort 0 | L8 | jac | V", 1 , 0, 8, "jac", 0.7, 8);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort0_l9_jac_v("IRBench cuda/cuda | q1 | sort 0 | L9 | jac | V", 1 , 0, 9, "jac", 0.7, 8);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort0_l10_jac_v("IRBench cuda/cuda | q1 | sort 0 | L10 | jac | V", 1 , 0, 10, "jac", 0.7, 8);

IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort1_l7_jac_v("IRBench cuda/cuda | q1 | sort 1 | L7 | jac | V", 1 , 1, 7, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort1_l8_jac_v("IRBench cuda/cuda | q1 | sort 1 | L8 | jac | V", 1 , 1, 8, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort1_l9_jac_v("IRBench cuda/cuda | q1 | sort 1 | L9 | jac | V", 1 , 1, 9, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort1_l10_jac_v("IRBench cuda/cuda | q1 | sort 1 | L10 | jac | V", 1 , 1, 10, "jac", 0.7, 15);

IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort2_l7_jac_v("IRBench cuda/cuda | q1 | sort 2 | L7 | jac | V", 1 , 2, 7, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort2_l8_jac_v("IRBench cuda/cuda | q1 | sort 2 | L8 | jac | V", 1 , 2, 8, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort2_l9_jac_v("IRBench cuda/cuda | q1 | sort 2 | L9 | jac | V", 1 , 2, 9, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort2_l10_jac_v("IRBench cuda/cuda | q1 | sort 2 | L10 | jac | V", 1 , 2, 10, "jac", 0.7, 15);

IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort4_l7_jac_v("IRBench cuda/cuda | q1 | sort 4 | L7 | jac | V", 1 , 4, 7, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort4_l8_jac_v("IRBench cuda/cuda | q1 | sort 4 | L8 | jac | V", 1 , 4, 8, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort4_l9_jac_v("IRBench cuda/cuda | q1 | sort 4 | L9 | jac | V", 1 , 4, 9, "jac", 0.7, 15);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q1_sort4_l10_jac_v("IRBench cuda/cuda | q1 | sort 4 | L10 | jac | V", 1 , 4, 10, "jac", 0.7, 15);

//q2
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q2_sort0_l6_jac_v("IRBench cuda/cuda | q2 | sort 0 | L6 | jac | V", 2 , 0, 6, "jac", 0.7, 30);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q2_sort0_l8_jac_v("IRBench cuda/cuda | q2 | sort 0 | L7 | jac | V", 2 , 0, 7, "jac", 0.7, 30);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q2_sort0_l9_jac_v("IRBench cuda/cuda | q2 | sort 0 | L8 | jac | V", 2 , 0, 8, "jac", 0.7, 30);
IRBench<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float, methods::CYCLE::V::STATIC, DenseVector<float> > cuda_q2_sort0_l10_jac_v("IRBench cuda/cuda | q2 | sort 0 | L9 | jac | V", 2 , 0, 9, "jac", 0.7, 30);

#endif
#endif
