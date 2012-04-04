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
#include <honei/math/mg.hh>
#include <honei/math/bicgstab.hh>
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


template <typename Tag_, typename CycleType_, typename PreconContType_>
class MGBench:
    public Benchmark
{
    private:
        unsigned long _element_type;
        unsigned long _sorting;
        unsigned long _levels;
        std::string _precon;
        double _damping;

    public:
        MGBench(const std::string & tag, unsigned long et, unsigned long s, unsigned long l, std::string p, double d) :
            Benchmark(tag)
        {
            register_tag(Tag_::name);
            _element_type = et;
            _sorting = s;
            _levels = l;
            _precon = p;
            _damping = d;
        }

        virtual void run()
        {
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson_advanced4/";

            if(_element_type == 2)
                file += "q2_";
            else if(_element_type == 3)
                file += "q1t_";

            file += "sort_";

            file += stringify(_sorting);

            file += "/";

            if(typeid(Tag_) == typeid(tags::CPU::MultiCore::SSE))
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


            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, PreconContType_, double>  data(MGUtil<Tag_,
                    SparseMatrixELL<double>,
                    DenseVector<double>,
                    SparseMatrixELL<double>,
                    PreconContType_,
                    io_formats::ELL,
                    io_formats::EXP,
                    double>::load_data(file, _levels, _damping, _precon));
            MGUtil<Tag_,
                SparseMatrixELL<double>,
                DenseVector<double>,
                SparseMatrixELL<double>,
                PreconContType_,
                io_formats::ELL,
                io_formats::EXP,
                double>::configure(data, 100, 150, 20, 20, 1, double(1e-8));

            OperatorList ol(
                    MGCycleCreation<Tag_,
                    CycleType_,
                    BiCGStabSolver<Tag_, methods::VAR>,
                    RISmoother<Tag_>,
                    Restriction<Tag_, methods::PROLMAT>,
                    Prolongation<Tag_, methods::PROLMAT>,
                    double>::value(data)
                    );

            BENCHMARK(
                    (MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data, ol));
#ifdef HONEI_CUDA
                    if (Tag_::tag_value == tags::tv_gpu_cuda)
                    cuda::GPUPool::instance()->flush();
#endif
                    );

            evaluate();

            std::cout << data.used_iters << " iterations used." << std::endl;
        }
};

#ifdef HONEI_SSE
//JAC
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_v("MGBench mcsse | q1 | sort 0 | L2 | jac | V", 1 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l8_jac_v("MGBench mcsse | q1 | sort 0 | L3 | jac | V", 1 , 0, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l9_jac_v("MGBench mcsse | q1 | sort 0 | L4 | jac | V", 1 , 0, 4, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_v("MGBench mcsse | q1 | sort 1 | L2 | jac | V", 1 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l8_jac_v("MGBench mcsse | q1 | sort 1 | L3 | jac | V", 1 , 1, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l9_jac_v("MGBench mcsse | q1 | sort 1 | L4 | jac | V", 1 , 1, 4, "jac", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_v("MGBench mcsse | q1 | sort 3 | L2 | jac | V", 1 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l8_jac_v("MGBench mcsse | q1 | sort 3 | L3 | jac | V", 1 , 3, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l9_jac_v("MGBench mcsse | q1 | sort 3 | L4 | jac | V", 1 , 3, 4, "jac", 0.5);






//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l7_jac_v("MGBench mcsse | q2 | sort 0 | L2 | jac | V", 2 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l8_jac_v("MGBench mcsse | q2 | sort 0 | L3 | jac | V", 2 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l7_jac_v("MGBench mcsse | q2 | sort 1 | L2 | jac | V", 2 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l8_jac_v("MGBench mcsse | q2 | sort 1 | L3 | jac | V", 2 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l7_jac_v("MGBench mcsse | q2 | sort 3 | L2 | jac | V", 2 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l8_jac_v("MGBench mcsse | q2 | sort 3 | L3 | jac | V", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l7_jac_v("MGBench mcsse | q1t | sort 0 | L2 | jac | V", 3 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l8_jac_v("MGBench mcsse | q1t | sort 0 | L3 | jac | V", 3 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l7_jac_v("MGBench mcsse | q1t | sort 1 | L2 | jac | V", 3 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l8_jac_v("MGBench mcsse | q1t | sort 1 | L3 | jac | V", 3 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l7_jac_v("MGBench mcsse | q1t | sort 3 | L2 | jac | V", 3 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l8_jac_v("MGBench mcsse | q1t | sort 3 | L3 | jac | V", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_grote_v("MGBench mcsse | q1 | sort 0 | L2 | spai_grote | V", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_spai_grote_v("MGBench mcsse | q1 | sort 0 | L3 | spai_grote | V", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_spai_grote_v("MGBench mcsse | q1 | sort 0 | L4 | spai_grote | V", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_grote_v("MGBench mcsse | q1 | sort 1 | L2 | spai_grote | V", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_spai_grote_v("MGBench mcsse | q1 | sort 1 | L3 | spai_grote | V", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_spai_grote_v("MGBench mcsse | q1 | sort 1 | L4 | spai_grote | V", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_grote_v("MGBench mcsse | q1 | sort 3 | L2 | spai_grote | V", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_spai_grote_v("MGBench mcsse | q1 | sort 3 | L3 | spai_grote | V", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_spai_grote_v("MGBench mcsse | q1 | sort 3 | L4 | spai_grote | V", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_spai_grote_v("MGBench mcsse | q2 | sort 0 | L2 | spai_grote | V", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_spai_grote_v("MGBench mcsse | q2 | sort 0 | L3 | spai_grote | V", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_spai_grote_v("MGBench mcsse | q2 | sort 1 | L2 | spai_grote | V", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_spai_grote_v("MGBench mcsse | q2 | sort 1 | L3 | spai_grote | V", 2 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_spai_grote_v("MGBench mcsse | q2 | sort 3 | L2 | spai_grote | V", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_spai_grote_v("MGBench mcsse | q2 | sort 3 | L3 | spai_grote | V", 2 , 3, 3, "spai_grote", 1.);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_spai_grote_v("MGBench mcsse | q1t | sort 0 | L2 | spai_grote | V", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_spai_grote_v("MGBench mcsse | q1t | sort 0 | L3 | spai_grote | V", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_spai_grote_v("MGBench mcsse | q1t | sort 1 | L2 | spai_grote | V", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_spai_grote_v("MGBench mcsse | q1t | sort 1 | L3 | spai_grote | V", 3 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_spai_grote_v("MGBench mcsse | q1t | sort 3 | L2 | spai_grote | V", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_spai_grote_v("MGBench mcsse | q1t | sort 3 | L3 | spai_grote | V", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_v("MGBench mcsse | q1 | sort 0 | L2 | sainv | V", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_sainv_v("MGBench mcsse | q1 | sort 0 | L3 | sainv | V", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_sainv_v("MGBench mcsse | q1 | sort 0 | L4 | sainv | V", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_v("MGBench mcsse | q1 | sort 1 | L2 | sainv | V", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_sainv_v("MGBench mcsse | q1 | sort 1 | L3 | sainv | V", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_sainv_v("MGBench mcsse | q1 | sort 1 | L4 | sainv | V", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_v("MGBench mcsse | q1 | sort 3 | L2 | sainv | V", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_sainv_v("MGBench mcsse | q1 | sort 3 | L3 | sainv | V", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_sainv_v("MGBench mcsse | q1 | sort 3 | L4 | sainv | V", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_sainv_v("MGBench mcsse | q2 | sort 0 | L2 | sainv | V", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_sainv_v("MGBench mcsse | q2 | sort 0 | L3 | sainv | V", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_sainv_v("MGBench mcsse | q2 | sort 1 | L2 | sainv | V", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_sainv_v("MGBench mcsse | q2 | sort 1 | L3 | sainv | V", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_sainv_v("MGBench mcsse | q2 | sort 3 | L2 | sainv | V", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_sainv_v("MGBench mcsse | q2 | sort 3 | L3 | sainv | V", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_sainv_v("MGBench mcsse | q1t | sort 0 | L2 | sainv | V", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_sainv_v("MGBench mcsse | q1t | sort 0 | L3 | sainv | V", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_sainv_v("MGBench mcsse | q1t | sort 1 | L2 | sainv | V", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_sainv_v("MGBench mcsse | q1t | sort 1 | L3 | sainv | V", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_sainv_v("MGBench mcsse | q1t | sort 3 | L2 | sainv | V", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_sainv_v("MGBench mcsse | q1t | sort 3 | L3 | sainv | V", 3 , 3, 3, "sainv", 0.5);




#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_v("MGBench cuda | q1 | sort 0 | L2 | jac | V", 1 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l8_jac_v("MGBench cuda | q1 | sort 0 | L3 | jac | V", 1 , 0, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l9_jac_v("MGBench cuda | q1 | sort 0 | L4 | jac | V", 1 , 0, 4, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_v("MGBench cuda | q1 | sort 1 | L2 | jac | V", 1 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l8_jac_v("MGBench cuda | q1 | sort 1 | L3 | jac | V", 1 , 1, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l9_jac_v("MGBench cuda | q1 | sort 1 | L4 | jac | V", 1 , 1, 4, "jac", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_v("MGBench cuda | q1 | sort 3 | L2 | jac | V", 1 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l8_jac_v("MGBench cuda | q1 | sort 3 | L3 | jac | V", 1 , 3, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l9_jac_v("MGBench cuda | q1 | sort 3 | L4 | jac | V", 1 , 3, 4, "jac", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort0_l7_jac_v("MGBench cuda | q2 | sort 0 | L2 | jac | V", 2 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort0_l8_jac_v("MGBench cuda | q2 | sort 0 | L3 | jac | V", 2 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort1_l7_jac_v("MGBench cuda | q2 | sort 1 | L2 | jac | V", 2 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort1_l8_jac_v("MGBench cuda | q2 | sort 1 | L3 | jac | V", 2 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort3_l7_jac_v("MGBench cuda | q2 | sort 3 | L2 | jac | V", 2 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort3_l8_jac_v("MGBench cuda | q2 | sort 3 | L3 | jac | V", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l7_jac_v("MGBench cuda | q1t | sort 0 | L2 | jac | V", 3 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l8_jac_v("MGBench cuda | q1t | sort 0 | L3 | jac | V", 3 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l7_jac_v("MGBench cuda | q1t | sort 1 | L2 | jac | V", 3 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l8_jac_v("MGBench cuda | q1t | sort 1 | L3 | jac | V", 3 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l7_jac_v("MGBench cuda | q1t | sort 3 | L2 | jac | V", 3 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l8_jac_v("MGBench cuda | q1t | sort 3 | L3 | jac | V", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_grote_v("MGBench cuda | q1 | sort 0 | L2 | spai_grote | V", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_spai_grote_v("MGBench cuda | q1 | sort 0 | L3 | spai_grote | V", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_spai_grote_v("MGBench cuda | q1 | sort 0 | L4 | spai_grote | V", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_grote_v("MGBench cuda | q1 | sort 1 | L2 | spai_grote | V", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_spai_grote_v("MGBench cuda | q1 | sort 1 | L3 | spai_grote | V", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_spai_grote_v("MGBench cuda | q1 | sort 1 | L4 | spai_grote | V", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_grote_v("MGBench cuda | q1 | sort 3 | L2 | spai_grote | V", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_spai_grote_v("MGBench cuda | q1 | sort 3 | L3 | spai_grote | V", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_spai_grote_v("MGBench cuda | q1 | sort 3 | L4 | spai_grote | V", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_spai_grote_v("MGBench cuda | q2 | sort 0 | L2 | spai_grote | V", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_spai_grote_v("MGBench cuda | q2 | sort 0 | L3 | spai_grote | V", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_spai_grote_v("MGBench cuda | q2 | sort 1 | L2 | spai_grote | V", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_spai_grote_v("MGBench cuda | q2 | sort 1 | L3 | spai_grote | V", 2 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_spai_grote_v("MGBench cuda | q2 | sort 3 | L2 | spai_grote | V", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_spai_grote_v("MGBench cuda | q2 | sort 3 | L3 | spai_grote | V", 2 , 3, 3, "spai_grote", 1.);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_spai_grote_v("MGBench cuda | q1t | sort 0 | L2 | spai_grote | V", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_spai_grote_v("MGBench cuda | q1t | sort 0 | L3 | spai_grote | V", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_spai_grote_v("MGBench cuda | q1t | sort 1 | L2 | spai_grote | V", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_spai_grote_v("MGBench cuda | q1t | sort 1 | L3 | spai_grote | V", 3 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_spai_grote_v("MGBench cuda | q1t | sort 3 | L2 | spai_grote | V", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_spai_grote_v("MGBench cuda | q1t | sort 3 | L3 | spai_grote | V", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_v("MGBench cuda | q1 | sort 0 | L2 | sainv | V", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_sainv_v("MGBench cuda | q1 | sort 0 | L3 | sainv | V", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_sainv_v("MGBench cuda | q1 | sort 0 | L4 | sainv | V", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_v("MGBench cuda | q1 | sort 1 | L2 | sainv | V", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_sainv_v("MGBench cuda | q1 | sort 1 | L3 | sainv | V", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_sainv_v("MGBench cuda | q1 | sort 1 | L4 | sainv | V", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_v("MGBench cuda | q1 | sort 3 | L2 | sainv | V", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_sainv_v("MGBench cuda | q1 | sort 3 | L3 | sainv | V", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_sainv_v("MGBench cuda | q1 | sort 3 | L4 | sainv | V", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_sainv_v("MGBench cuda | q2 | sort 0 | L2 | sainv | V", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_sainv_v("MGBench cuda | q2 | sort 0 | L3 | sainv | V", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_sainv_v("MGBench cuda | q2 | sort 1 | L2 | sainv | V", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_sainv_v("MGBench cuda | q2 | sort 1 | L3 | sainv | V", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_sainv_v("MGBench cuda | q2 | sort 3 | L2 | sainv | V", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_sainv_v("MGBench cuda | q2 | sort 3 | L3 | sainv | V", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_sainv_v("MGBench cuda | q1t | sort 0 | L2 | sainv | V", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_sainv_v("MGBench cuda | q1t | sort 0 | L3 | sainv | V", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_sainv_v("MGBench cuda | q1t | sort 1 | L2 | sainv | V", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_sainv_v("MGBench cuda | q1t | sort 1 | L3 | sainv | V", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_sainv_v("MGBench cuda | q1t | sort 3 | L2 | sainv | V", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_sainv_v("MGBench cuda | q1t | sort 3 | L3 | sainv | V", 3 , 3, 3, "sainv", 0.5);





#endif

//W-cycle
#ifdef HONEI_SSE
//JAC
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_w("MGBench mcsse | q1 | sort 0 | L2 | jac | W", 1 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l8_jac_w("MGBench mcsse | q1 | sort 0 | L3 | jac | W", 1 , 0, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l9_jac_w("MGBench mcsse | q1 | sort 0 | L4 | jac | W", 1 , 0, 4, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_w("MGBench mcsse | q1 | sort 1 | L2 | jac | W", 1 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l8_jac_w("MGBench mcsse | q1 | sort 1 | L3 | jac | W", 1 , 1, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l9_jac_w("MGBench mcsse | q1 | sort 1 | L4 | jac | W", 1 , 1, 4, "jac", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_w("MGBench mcsse | q1 | sort 3 | L2 | jac | W", 1 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l8_jac_w("MGBench mcsse | q1 | sort 3 | L3 | jac | W", 1 , 3, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l9_jac_w("MGBench mcsse | q1 | sort 3 | L4 | jac | W", 1 , 3, 4, "jac", 0.5);






//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort0_l7_jac_w("MGBench mcsse | q2 | sort 0 | L2 | jac | W", 2 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort0_l8_jac_w("MGBench mcsse | q2 | sort 0 | L3 | jac | W", 2 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort1_l7_jac_w("MGBench mcsse | q2 | sort 1 | L2 | jac | W", 2 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort1_l8_jac_w("MGBench mcsse | q2 | sort 1 | L3 | jac | W", 2 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort3_l7_jac_w("MGBench mcsse | q2 | sort 3 | L2 | jac | W", 2 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort3_l8_jac_w("MGBench mcsse | q2 | sort 3 | L3 | jac | W", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort0_l7_jac_w("MGBench mcsse | q1t | sort 0 | L2 | jac | W", 3 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort0_l8_jac_w("MGBench mcsse | q1t | sort 0 | L3 | jac | W", 3 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort1_l7_jac_w("MGBench mcsse | q1t | sort 1 | L2 | jac | W", 3 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort1_l8_jac_w("MGBench mcsse | q1t | sort 1 | L3 | jac | W", 3 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort3_l7_jac_w("MGBench mcsse | q1t | sort 3 | L2 | jac | W", 3 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort3_l8_jac_w("MGBench mcsse | q1t | sort 3 | L3 | jac | W", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_grote_w("MGBench mcsse | q1 | sort 0 | L2 | spai_grote | W", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_spai_grote_w("MGBench mcsse | q1 | sort 0 | L3 | spai_grote | W", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_spai_grote_w("MGBench mcsse | q1 | sort 0 | L4 | spai_grote | W", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_grote_w("MGBench mcsse | q1 | sort 1 | L2 | spai_grote | W", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_spai_grote_w("MGBench mcsse | q1 | sort 1 | L3 | spai_grote | W", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_spai_grote_w("MGBench mcsse | q1 | sort 1 | L4 | spai_grote | W", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_grote_w("MGBench mcsse | q1 | sort 3 | L2 | spai_grote | W", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_spai_grote_w("MGBench mcsse | q1 | sort 3 | L3 | spai_grote | W", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_spai_grote_w("MGBench mcsse | q1 | sort 3 | L4 | spai_grote | W", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_spai_grote_w("MGBench mcsse | q2 | sort 0 | L2 | spai_grote | W", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_spai_grote_w("MGBench mcsse | q2 | sort 0 | L3 | spai_grote | W", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_spai_grote_w("MGBench mcsse | q2 | sort 1 | L2 | spai_grote | W", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_spai_grote_w("MGBench mcsse | q2 | sort 1 | L3 | spai_grote | W", 2 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_spai_grote_w("MGBench mcsse | q2 | sort 3 | L2 | spai_grote | W", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_spai_grote_w("MGBench mcsse | q2 | sort 3 | L3 | spai_grote | W", 2 , 3, 3, "spai_grote", 1.);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_spai_grote_w("MGBench mcsse | q1t | sort 0 | L2 | spai_grote | W", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_spai_grote_w("MGBench mcsse | q1t | sort 0 | L3 | spai_grote | W", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_spai_grote_w("MGBench mcsse | q1t | sort 1 | L2 | spai_grote | W", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_spai_grote_w("MGBench mcsse | q1t | sort 1 | L3 | spai_grote | W", 3 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_spai_grote_w("MGBench mcsse | q1t | sort 3 | L2 | spai_grote | W", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_spai_grote_w("MGBench mcsse | q1t | sort 3 | L3 | spai_grote | W", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_w("MGBench mcsse | q1 | sort 0 | L2 | sainv | W", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_sainv_w("MGBench mcsse | q1 | sort 0 | L3 | sainv | W", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_sainv_w("MGBench mcsse | q1 | sort 0 | L4 | sainv | W", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_w("MGBench mcsse | q1 | sort 1 | L2 | sainv | W", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_sainv_w("MGBench mcsse | q1 | sort 1 | L3 | sainv | W", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_sainv_w("MGBench mcsse | q1 | sort 1 | L4 | sainv | W", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_w("MGBench mcsse | q1 | sort 3 | L2 | sainv | W", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_sainv_w("MGBench mcsse | q1 | sort 3 | L3 | sainv | W", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_sainv_w("MGBench mcsse | q1 | sort 3 | L4 | sainv | W", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_sainv_w("MGBench mcsse | q2 | sort 0 | L2 | sainv | W", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_sainv_w("MGBench mcsse | q2 | sort 0 | L3 | sainv | W", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_sainv_w("MGBench mcsse | q2 | sort 1 | L2 | sainv | W", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_sainv_w("MGBench mcsse | q2 | sort 1 | L3 | sainv | W", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_sainv_w("MGBench mcsse | q2 | sort 3 | L2 | sainv | W", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_sainv_w("MGBench mcsse | q2 | sort 3 | L3 | sainv | W", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_sainv_w("MGBench mcsse | q1t | sort 0 | L2 | sainv | W", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_sainv_w("MGBench mcsse | q1t | sort 0 | L3 | sainv | W", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_sainv_w("MGBench mcsse | q1t | sort 1 | L2 | sainv | W", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_sainv_w("MGBench mcsse | q1t | sort 1 | L3 | sainv | W", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_sainv_w("MGBench mcsse | q1t | sort 3 | L2 | sainv | W", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_sainv_w("MGBench mcsse | q1t | sort 3 | L3 | sainv | W", 3 , 3, 3, "sainv", 0.5);




#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_w("MGBench cuda | q1 | sort 0 | L2 | jac | W", 1 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l8_jac_w("MGBench cuda | q1 | sort 0 | L3 | jac | W", 1 , 0, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l9_jac_w("MGBench cuda | q1 | sort 0 | L4 | jac | W", 1 , 0, 4, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_w("MGBench cuda | q1 | sort 1 | L2 | jac | W", 1 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l8_jac_w("MGBench cuda | q1 | sort 1 | L3 | jac | W", 1 , 1, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l9_jac_w("MGBench cuda | q1 | sort 1 | L4 | jac | W", 1 , 1, 4, "jac", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_w("MGBench cuda | q1 | sort 3 | L2 | jac | W", 1 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l8_jac_w("MGBench cuda | q1 | sort 3 | L3 | jac | W", 1 , 3, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l9_jac_w("MGBench cuda | q1 | sort 3 | L4 | jac | W", 1 , 3, 4, "jac", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort0_l7_jac_w("MGBench cuda | q2 | sort 0 | L2 | jac | W", 2 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort0_l8_jac_w("MGBench cuda | q2 | sort 0 | L3 | jac | W", 2 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort1_l7_jac_w("MGBench cuda | q2 | sort 1 | L2 | jac | W", 2 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort1_l8_jac_w("MGBench cuda | q2 | sort 1 | L3 | jac | W", 2 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort3_l7_jac_w("MGBench cuda | q2 | sort 3 | L2 | jac | W", 2 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort3_l8_jac_w("MGBench cuda | q2 | sort 3 | L3 | jac | W", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort0_l7_jac_w("MGBench cuda | q1t | sort 0 | L2 | jac | W", 3 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort0_l8_jac_w("MGBench cuda | q1t | sort 0 | L3 | jac | W", 3 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort1_l7_jac_w("MGBench cuda | q1t | sort 1 | L2 | jac | W", 3 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort1_l8_jac_w("MGBench cuda | q1t | sort 1 | L3 | jac | W", 3 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort3_l7_jac_w("MGBench cuda | q1t | sort 3 | L2 | jac | W", 3 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort3_l8_jac_w("MGBench cuda | q1t | sort 3 | L3 | jac | W", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_grote_w("MGBench cuda | q1 | sort 0 | L2 | spai_grote | W", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_spai_grote_w("MGBench cuda | q1 | sort 0 | L3 | spai_grote | W", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_spai_grote_w("MGBench cuda | q1 | sort 0 | L4 | spai_grote | W", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_grote_w("MGBench cuda | q1 | sort 1 | L2 | spai_grote | W", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_spai_grote_w("MGBench cuda | q1 | sort 1 | L3 | spai_grote | W", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_spai_grote_w("MGBench cuda | q1 | sort 1 | L4 | spai_grote | W", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_grote_w("MGBench cuda | q1 | sort 3 | L2 | spai_grote | W", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_spai_grote_w("MGBench cuda | q1 | sort 3 | L3 | spai_grote | W", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_spai_grote_w("MGBench cuda | q1 | sort 3 | L4 | spai_grote | W", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_spai_grote_w("MGBench cuda | q2 | sort 0 | L2 | spai_grote | W", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_spai_grote_w("MGBench cuda | q2 | sort 0 | L3 | spai_grote | W", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_spai_grote_w("MGBench cuda | q2 | sort 1 | L2 | spai_grote | W", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_spai_grote_w("MGBench cuda | q2 | sort 1 | L3 | spai_grote | W", 2 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_spai_grote_w("MGBench cuda | q2 | sort 3 | L2 | spai_grote | W", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_spai_grote_w("MGBench cuda | q2 | sort 3 | L3 | spai_grote | W", 2 , 3, 3, "spai_grote", 1.);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_spai_grote_w("MGBench cuda | q1t | sort 0 | L2 | spai_grote | W", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_spai_grote_w("MGBench cuda | q1t | sort 0 | L3 | spai_grote | W", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_spai_grote_w("MGBench cuda | q1t | sort 1 | L2 | spai_grote | W", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_spai_grote_w("MGBench cuda | q1t | sort 1 | L3 | spai_grote | W", 3 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_spai_grote_w("MGBench cuda | q1t | sort 3 | L2 | spai_grote | W", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_spai_grote_w("MGBench cuda | q1t | sort 3 | L3 | spai_grote | W", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_w("MGBench cuda | q1 | sort 0 | L2 | sainv | W", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_sainv_w("MGBench cuda | q1 | sort 0 | L3 | sainv | W", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_sainv_w("MGBench cuda | q1 | sort 0 | L4 | sainv | W", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_w("MGBench cuda | q1 | sort 1 | L2 | sainv | W", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_sainv_w("MGBench cuda | q1 | sort 1 | L3 | sainv | W", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_sainv_w("MGBench cuda | q1 | sort 1 | L4 | sainv | W", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_w("MGBench cuda | q1 | sort 3 | L2 | sainv | W", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_sainv_w("MGBench cuda | q1 | sort 3 | L3 | sainv | W", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_sainv_w("MGBench cuda | q1 | sort 3 | L4 | sainv | W", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_sainv_w("MGBench cuda | q2 | sort 0 | L2 | sainv | W", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_sainv_w("MGBench cuda | q2 | sort 0 | L3 | sainv | W", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_sainv_w("MGBench cuda | q2 | sort 1 | L2 | sainv | W", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_sainv_w("MGBench cuda | q2 | sort 1 | L3 | sainv | W", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_sainv_w("MGBench cuda | q2 | sort 3 | L2 | sainv | W", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_sainv_w("MGBench cuda | q2 | sort 3 | L3 | sainv | W", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_sainv_w("MGBench cuda | q1t | sort 0 | L2 | sainv | W", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_sainv_w("MGBench cuda | q1t | sort 0 | L3 | sainv | W", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_sainv_w("MGBench cuda | q1t | sort 1 | L2 | sainv | W", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_sainv_w("MGBench cuda | q1t | sort 1 | L3 | sainv | W", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_sainv_w("MGBench cuda | q1t | sort 3 | L2 | sainv | W", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_sainv_w("MGBench cuda | q1t | sort 3 | L3 | sainv | W", 3 , 3, 3, "sainv", 0.5);





#endif


//F-cycle
#ifdef HONEI_SSE
//JAC
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_f("MGBench mcsse | q1 | sort 0 | L2 | jac | F", 1 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l8_jac_f("MGBench mcsse | q1 | sort 0 | L3 | jac | F", 1 , 0, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l9_jac_f("MGBench mcsse | q1 | sort 0 | L4 | jac | F", 1 , 0, 4, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_f("MGBench mcsse | q1 | sort 1 | L2 | jac | F", 1 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l8_jac_f("MGBench mcsse | q1 | sort 1 | L3 | jac | F", 1 , 1, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l9_jac_f("MGBench mcsse | q1 | sort 1 | L4 | jac | F", 1 , 1, 4, "jac", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_f("MGBench mcsse | q1 | sort 3 | L2 | jac | F", 1 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l8_jac_f("MGBench mcsse | q1 | sort 3 | L3 | jac | F", 1 , 3, 3, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l9_jac_f("MGBench mcsse | q1 | sort 3 | L4 | jac | F", 1 , 3, 4, "jac", 0.5);






//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l7_jac_f("MGBench mcsse | q2 | sort 0 | L2 | jac | F", 2 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l8_jac_f("MGBench mcsse | q2 | sort 0 | L3 | jac | F", 2 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l7_jac_f("MGBench mcsse | q2 | sort 1 | L2 | jac | F", 2 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l8_jac_f("MGBench mcsse | q2 | sort 1 | L3 | jac | F", 2 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l7_jac_f("MGBench mcsse | q2 | sort 3 | L2 | jac | F", 2 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l8_jac_f("MGBench mcsse | q2 | sort 3 | L3 | jac | F", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l7_jac_f("MGBench mcsse | q1t | sort 0 | L2 | jac | F", 3 , 0, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l8_jac_f("MGBench mcsse | q1t | sort 0 | L3 | jac | F", 3 , 0, 3, "jac", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l7_jac_f("MGBench mcsse | q1t | sort 1 | L2 | jac | F", 3 , 1, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l8_jac_f("MGBench mcsse | q1t | sort 1 | L3 | jac | F", 3 , 1, 3, "jac", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l7_jac_f("MGBench mcsse | q1t | sort 3 | L2 | jac | F", 3 , 3, 2, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l8_jac_f("MGBench mcsse | q1t | sort 3 | L3 | jac | F", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_grote_f("MGBench mcsse | q1 | sort 0 | L2 | spai_grote | F", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_spai_grote_f("MGBench mcsse | q1 | sort 0 | L3 | spai_grote | F", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_spai_grote_f("MGBench mcsse | q1 | sort 0 | L4 | spai_grote | F", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_grote_f("MGBench mcsse | q1 | sort 1 | L2 | spai_grote | F", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_spai_grote_f("MGBench mcsse | q1 | sort 1 | L3 | spai_grote | F", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_spai_grote_f("MGBench mcsse | q1 | sort 1 | L4 | spai_grote | F", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_grote_f("MGBench mcsse | q1 | sort 3 | L2 | spai_grote | F", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_spai_grote_f("MGBench mcsse | q1 | sort 3 | L3 | spai_grote | F", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_spai_grote_f("MGBench mcsse | q1 | sort 3 | L4 | spai_grote | F", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_spai_grote_f("MGBench mcsse | q2 | sort 0 | L2 | spai_grote | F", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_spai_grote_f("MGBench mcsse | q2 | sort 0 | L3 | spai_grote | F", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_spai_grote_f("MGBench mcsse | q2 | sort 1 | L2 | spai_grote | F", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_spai_grote_f("MGBench mcsse | q2 | sort 1 | L3 | spai_grote | F", 2 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_spai_grote_f("MGBench mcsse | q2 | sort 3 | L2 | spai_grote | F", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_spai_grote_f("MGBench mcsse | q2 | sort 3 | L3 | spai_grote | F", 2 , 3, 3, "spai_grote", 1.);





//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_spai_grote_f("MGBench mcsse | q1t | sort 0 | L2 | spai_grote | F", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_spai_grote_f("MGBench mcsse | q1t | sort 0 | L3 | spai_grote | F", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_spai_grote_f("MGBench mcsse | q1t | sort 1 | L2 | spai_grote | F", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_spai_grote_f("MGBench mcsse | q1t | sort 1 | L3 | spai_grote | F", 3 , 1, 3, "spai_grote", 1.);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_spai_grote_f("MGBench mcsse | q1t | sort 3 | L2 | spai_grote | F", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_spai_grote_f("MGBench mcsse | q1t | sort 3 | L3 | spai_grote | F", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_f("MGBench mcsse | q1 | sort 0 | L2 | sainv | F", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l8_sainv_f("MGBench mcsse | q1 | sort 0 | L3 | sainv | F", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l9_sainv_f("MGBench mcsse | q1 | sort 0 | L4 | sainv | F", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_f("MGBench mcsse | q1 | sort 1 | L2 | sainv | F", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l8_sainv_f("MGBench mcsse | q1 | sort 1 | L3 | sainv | F", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l9_sainv_f("MGBench mcsse | q1 | sort 1 | L4 | sainv | F", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_f("MGBench mcsse | q1 | sort 3 | L2 | sainv | F", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l8_sainv_f("MGBench mcsse | q1 | sort 3 | L3 | sainv | F", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l9_sainv_f("MGBench mcsse | q1 | sort 3 | L4 | sainv | F", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l7_sainv_f("MGBench mcsse | q2 | sort 0 | L2 | sainv | F", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l8_sainv_f("MGBench mcsse | q2 | sort 0 | L3 | sainv | F", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l7_sainv_f("MGBench mcsse | q2 | sort 1 | L2 | sainv | F", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l8_sainv_f("MGBench mcsse | q2 | sort 1 | L3 | sainv | F", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l7_sainv_f("MGBench mcsse | q2 | sort 3 | L2 | sainv | F", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l8_sainv_f("MGBench mcsse | q2 | sort 3 | L3 | sainv | F", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l7_sainv_f("MGBench mcsse | q1t | sort 0 | L2 | sainv | F", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l8_sainv_f("MGBench mcsse | q1t | sort 0 | L3 | sainv | F", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l7_sainv_f("MGBench mcsse | q1t | sort 1 | L2 | sainv | F", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l8_sainv_f("MGBench mcsse | q1t | sort 1 | L3 | sainv | F", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l7_sainv_f("MGBench mcsse | q1t | sort 3 | L2 | sainv | F", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l8_sainv_f("MGBench mcsse | q1t | sort 3 | L3 | sainv | F", 3 , 3, 3, "sainv", 0.5);




#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_f("MGBench cuda | q1 | sort 0 | L2 | jac | F", 1 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l8_jac_f("MGBench cuda | q1 | sort 0 | L3 | jac | F", 1 , 0, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l9_jac_f("MGBench cuda | q1 | sort 0 | L4 | jac | F", 1 , 0, 4, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_f("MGBench cuda | q1 | sort 1 | L2 | jac | F", 1 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l8_jac_f("MGBench cuda | q1 | sort 1 | L3 | jac | F", 1 , 1, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l9_jac_f("MGBench cuda | q1 | sort 1 | L4 | jac | F", 1 , 1, 4, "jac", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_f("MGBench cuda | q1 | sort 3 | L2 | jac | F", 1 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l8_jac_f("MGBench cuda | q1 | sort 3 | L3 | jac | F", 1 , 3, 3, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l9_jac_f("MGBench cuda | q1 | sort 3 | L4 | jac | F", 1 , 3, 4, "jac", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort0_l7_jac_f("MGBench cuda | q2 | sort 0 | L2 | jac | F", 2 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort0_l8_jac_f("MGBench cuda | q2 | sort 0 | L3 | jac | F", 2 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort1_l7_jac_f("MGBench cuda | q2 | sort 1 | L2 | jac | F", 2 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort1_l8_jac_f("MGBench cuda | q2 | sort 1 | L3 | jac | F", 2 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort3_l7_jac_f("MGBench cuda | q2 | sort 3 | L2 | jac | F", 2 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort3_l8_jac_f("MGBench cuda | q2 | sort 3 | L3 | jac | F", 2 , 3, 3, "jac", 0.5);





//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l7_jac_f("MGBench cuda | q1t | sort 0 | L2 | jac | F", 3 , 0, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l8_jac_f("MGBench cuda | q1t | sort 0 | L3 | jac | F", 3 , 0, 3, "jac", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l7_jac_f("MGBench cuda | q1t | sort 1 | L2 | jac | F", 3 , 1, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l8_jac_f("MGBench cuda | q1t | sort 1 | L3 | jac | F", 3 , 1, 3, "jac", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l7_jac_f("MGBench cuda | q1t | sort 3 | L2 | jac | F", 3 , 3, 2, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l8_jac_f("MGBench cuda | q1t | sort 3 | L3 | jac | F", 3 , 3, 3, "jac", 0.5);





//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_grote_f("MGBench cuda | q1 | sort 0 | L2 | spai_grote | F", 1 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_spai_grote_f("MGBench cuda | q1 | sort 0 | L3 | spai_grote | F", 1 , 0, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_spai_grote_f("MGBench cuda | q1 | sort 0 | L4 | spai_grote | F", 1 , 0, 4, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_grote_f("MGBench cuda | q1 | sort 1 | L2 | spai_grote | F", 1 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_spai_grote_f("MGBench cuda | q1 | sort 1 | L3 | spai_grote | F", 1 , 1, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_spai_grote_f("MGBench cuda | q1 | sort 1 | L4 | spai_grote | F", 1 , 1, 4, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_grote_f("MGBench cuda | q1 | sort 3 | L2 | spai_grote | F", 1 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_spai_grote_f("MGBench cuda | q1 | sort 3 | L3 | spai_grote | F", 1 , 3, 3, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_spai_grote_f("MGBench cuda | q1 | sort 3 | L4 | spai_grote | F", 1 , 3, 4, "spai_grote", 1.);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_spai_grote_f("MGBench cuda | q2 | sort 0 | L2 | spai_grote | F", 2 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_spai_grote_f("MGBench cuda | q2 | sort 0 | L3 | spai_grote | F", 2 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_spai_grote_f("MGBench cuda | q2 | sort 1 | L2 | spai_grote | F", 2 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_spai_grote_f("MGBench cuda | q2 | sort 1 | L3 | spai_grote | F", 2 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_spai_grote_f("MGBench cuda | q2 | sort 3 | L2 | spai_grote | F", 2 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_spai_grote_f("MGBench cuda | q2 | sort 3 | L3 | spai_grote | F", 2 , 3, 3, "spai_grote", 1.);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_spai_grote_f("MGBench cuda | q1t | sort 0 | L2 | spai_grote | F", 3 , 0, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_spai_grote_f("MGBench cuda | q1t | sort 0 | L3 | spai_grote | F", 3 , 0, 3, "spai_grote", 1.);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_spai_grote_f("MGBench cuda | q1t | sort 1 | L2 | spai_grote | F", 3 , 1, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_spai_grote_f("MGBench cuda | q1t | sort 1 | L3 | spai_grote | F", 3 , 1, 3, "spai_grote", 1.);






MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_spai_grote_f("MGBench cuda | q1t | sort 3 | L2 | spai_grote | F", 3 , 3, 2, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_spai_grote_f("MGBench cuda | q1t | sort 3 | L3 | spai_grote | F", 3 , 3, 3, "spai_grote", 1.);





//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_f("MGBench cuda | q1 | sort 0 | L2 | sainv | F", 1 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l8_sainv_f("MGBench cuda | q1 | sort 0 | L3 | sainv | F", 1 , 0, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l9_sainv_f("MGBench cuda | q1 | sort 0 | L4 | sainv | F", 1 , 0, 4, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_f("MGBench cuda | q1 | sort 1 | L2 | sainv | F", 1 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l8_sainv_f("MGBench cuda | q1 | sort 1 | L3 | sainv | F", 1 , 1, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l9_sainv_f("MGBench cuda | q1 | sort 1 | L4 | sainv | F", 1 , 1, 4, "sainv", 0.5);






MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_f("MGBench cuda | q1 | sort 3 | L2 | sainv | F", 1 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l8_sainv_f("MGBench cuda | q1 | sort 3 | L3 | sainv | F", 1 , 3, 3, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l9_sainv_f("MGBench cuda | q1 | sort 3 | L4 | sainv | F", 1 , 3, 4, "sainv", 0.5);





//q2

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l7_sainv_f("MGBench cuda | q2 | sort 0 | L2 | sainv | F", 2 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l8_sainv_f("MGBench cuda | q2 | sort 0 | L3 | sainv | F", 2 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l7_sainv_f("MGBench cuda | q2 | sort 1 | L2 | sainv | F", 2 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l8_sainv_f("MGBench cuda | q2 | sort 1 | L3 | sainv | F", 2 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l7_sainv_f("MGBench cuda | q2 | sort 3 | L2 | sainv | F", 2 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l8_sainv_f("MGBench cuda | q2 | sort 3 | L3 | sainv | F", 2 , 3, 3, "sainv", 0.5);






//q1t

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l7_sainv_f("MGBench cuda | q1t | sort 0 | L2 | sainv | F", 3 , 0, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l8_sainv_f("MGBench cuda | q1t | sort 0 | L3 | sainv | F", 3 , 0, 3, "sainv", 0.5);



MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l7_sainv_f("MGBench cuda | q1t | sort 1 | L2 | sainv | F", 3 , 1, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l8_sainv_f("MGBench cuda | q1t | sort 1 | L3 | sainv | F", 3 , 1, 3, "sainv", 0.5);







MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l7_sainv_f("MGBench cuda | q1t | sort 3 | L2 | sainv | F", 3 , 3, 2, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l8_sainv_f("MGBench cuda | q1t | sort 3 | L3 | sainv | F", 3 , 3, 3, "sainv", 0.5);





#endif
