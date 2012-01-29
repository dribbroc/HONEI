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
#include <honei/math/superlu.hh>

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


            ///SET ADAPTIVE SMOOTHING HERE (Krylov Smoother)
            bool adaptive = true;

            _damping = adaptive ? double(1) : _damping;

            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, PreconContType_, double >  data(MGUtil<Tag_,
                    SparseMatrixELL<double>,
                    DenseVector<double>,
                    SparseMatrixELL<double>,
                    PreconContType_,
                    io_formats::ELL,
                    io_formats::EXP,
                    double,
                    BiCGStabSmoother<Tag_> >::load_data(file, _levels, _damping, _precon));
            MGUtil<Tag_,
                SparseMatrixELL<double>,
                DenseVector<double>,
                SparseMatrixELL<double>,
                PreconContType_,
                io_formats::ELL,
                io_formats::EXP,
                double,
                BiCGStabSmoother<Tag_> >::configure(data, 100, 10, 2, 2, 1, double(1e-8));

            OperatorList ol(
                    MGCycleCreation<Tag_,
                    CycleType_,
                    //SuperLU,
                    BiCGStab<Tag_, methods::VAR>,
                    //CG<Tag_, methods::NONE>,
                    //RISmoother<Tag_>,
                    BiCGStabSmoother<Tag_>,
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
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l5_jac_v("MGBench mcsse | q1 | sort 0 | L5 | jac | V", 1 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l6_jac_v("MGBench mcsse | q1 | sort 0 | L6 | jac | V", 1 , 0, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_v("MGBench mcsse | q1 | sort 0 | L7 | jac | V", 1 , 0, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l5_jac_v("MGBench mcsse | q1 | sort 1 | L5 | jac | V", 1 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l6_jac_v("MGBench mcsse | q1 | sort 1 | L6 | jac | V", 1 , 1, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_v("MGBench mcsse | q1 | sort 1 | L7 | jac | V", 1 , 1, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l5_jac_v("MGBench mcsse | q1 | sort 2 | L5 | jac | V", 1 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l6_jac_v("MGBench mcsse | q1 | sort 2 | L6 | jac | V", 1 , 2, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l7_jac_v("MGBench mcsse | q1 | sort 2 | L7 | jac | V", 1 , 2, 7, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l5_jac_v("MGBench mcsse | q1 | sort 3 | L5 | jac | V", 1 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l6_jac_v("MGBench mcsse | q1 | sort 3 | L6 | jac | V", 1 , 3, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_v("MGBench mcsse | q1 | sort 3 | L7 | jac | V", 1 , 3, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l5_jac_v("MGBench mcsse | q1 | sort 4 | L5 | jac | V", 1 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l6_jac_v("MGBench mcsse | q1 | sort 4 | L6 | jac | V", 1 , 4, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l7_jac_v("MGBench mcsse | q1 | sort 4 | L7 | jac | V", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l4_jac_v("MGBench mcsse | q2 | sort 0 | L4 | jac | V", 2 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l5_jac_v("MGBench mcsse | q2 | sort 0 | L5 | jac | V", 2 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l6_jac_v("MGBench mcsse | q2 | sort 0 | L6 | jac | V", 2 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l4_jac_v("MGBench mcsse | q2 | sort 1 | L4 | jac | V", 2 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l5_jac_v("MGBench mcsse | q2 | sort 1 | L5 | jac | V", 2 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l6_jac_v("MGBench mcsse | q2 | sort 1 | L6 | jac | V", 2 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l4_jac_v("MGBench mcsse | q2 | sort 2 | L4 | jac | V", 2 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l5_jac_v("MGBench mcsse | q2 | sort 2 | L5 | jac | V", 2 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l6_jac_v("MGBench mcsse | q2 | sort 2 | L6 | jac | V", 2 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l4_jac_v("MGBench mcsse | q2 | sort 3 | L4 | jac | V", 2 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l5_jac_v("MGBench mcsse | q2 | sort 3 | L5 | jac | V", 2 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l6_jac_v("MGBench mcsse | q2 | sort 3 | L6 | jac | V", 2 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l4_jac_v("MGBench mcsse | q2 | sort 4 | L4 | jac | V", 2 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l5_jac_v("MGBench mcsse | q2 | sort 4 | L5 | jac | V", 2 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l6_jac_v("MGBench mcsse | q2 | sort 4 | L6 | jac | V", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l4_jac_v("MGBench mcsse | q1t | sort 0 | L4 | jac | V", 3 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l5_jac_v("MGBench mcsse | q1t | sort 0 | L5 | jac | V", 3 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l6_jac_v("MGBench mcsse | q1t | sort 0 | L6 | jac | V", 3 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l4_jac_v("MGBench mcsse | q1t | sort 1 | L4 | jac | V", 3 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l5_jac_v("MGBench mcsse | q1t | sort 1 | L5 | jac | V", 3 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l6_jac_v("MGBench mcsse | q1t | sort 1 | L6 | jac | V", 3 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l4_jac_v("MGBench mcsse | q1t | sort 2 | L4 | jac | V", 3 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l5_jac_v("MGBench mcsse | q1t | sort 2 | L5 | jac | V", 3 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l6_jac_v("MGBench mcsse | q1t | sort 2 | L6 | jac | V", 3 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l4_jac_v("MGBench mcsse | q1t | sort 3 | L4 | jac | V", 3 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l5_jac_v("MGBench mcsse | q1t | sort 3 | L5 | jac | V", 3 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l6_jac_v("MGBench mcsse | q1t | sort 3 | L6 | jac | V", 3 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l4_jac_v("MGBench mcsse | q1t | sort 4 | L4 | jac | V", 3 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l5_jac_v("MGBench mcsse | q1t | sort 4 | L5 | jac | V", 3 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l6_jac_v("MGBench mcsse | q1t | sort 4 | L6 | jac | V", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_spai_v("MGBench mcsse | q1 | sort 0 | L5 | spai | V", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_spai_v("MGBench mcsse | q1 | sort 0 | L6 | spai | V", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_v("MGBench mcsse | q1 | sort 0 | L7 | spai | V", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_spai_v("MGBench mcsse | q1 | sort 1 | L5 | spai | V", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_spai_v("MGBench mcsse | q1 | sort 1 | L6 | spai | V", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_v("MGBench mcsse | q1 | sort 1 | L7 | spai | V", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_spai_v("MGBench mcsse | q1 | sort 2 | L5 | spai | V", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_spai_v("MGBench mcsse | q1 | sort 2 | L6 | spai | V", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_spai_v("MGBench mcsse | q1 | sort 2 | L7 | spai | V", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_spai_v("MGBench mcsse | q1 | sort 3 | L5 | spai | V", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_spai_v("MGBench mcsse | q1 | sort 3 | L6 | spai | V", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_v("MGBench mcsse | q1 | sort 3 | L7 | spai | V", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_spai_v("MGBench mcsse | q1 | sort 4 | L5 | spai | V", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_spai_v("MGBench mcsse | q1 | sort 4 | L6 | spai | V", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_spai_v("MGBench mcsse | q1 | sort 4 | L7 | spai | V", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_spai_v("MGBench mcsse | q2 | sort 0 | L4 | spai | V", 2 , 0, 4, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_spai_v("MGBench mcsse | q2 | sort 0 | L5 | spai | V", 2 , 0, 5, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_spai_v("MGBench mcsse | q2 | sort 0 | L6 | spai | V", 2 , 0, 6, "spai_grote", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_spai_v("MGBench mcsse | q2 | sort 1 | L4 | spai | V", 2 , 1, 4, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_spai_v("MGBench mcsse | q2 | sort 1 | L5 | spai | V", 2 , 1, 5, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_spai_v("MGBench mcsse | q2 | sort 1 | L6 | spai | V", 2 , 1, 6, "spai_grote", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_spai_v("MGBench mcsse | q2 | sort 2 | L4 | spai | V", 2 , 2, 4, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_spai_v("MGBench mcsse | q2 | sort 2 | L5 | spai | V", 2 , 2, 5, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_spai_v("MGBench mcsse | q2 | sort 2 | L6 | spai | V", 2 , 2, 6, "spai_grote", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_spai_v("MGBench mcsse | q2 | sort 3 | L4 | spai | V", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_spai_v("MGBench mcsse | q2 | sort 3 | L5 | spai | V", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_spai_v("MGBench mcsse | q2 | sort 3 | L6 | spai | V", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_spai_v("MGBench mcsse | q2 | sort 4 | L4 | spai | V", 2 , 4, 4, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_spai_v("MGBench mcsse | q2 | sort 4 | L5 | spai | V", 2 , 4, 5, "spai_grote", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_spai_v("MGBench mcsse | q2 | sort 4 | L6 | spai | V", 2 , 4, 6, "spai_grote", 0.5);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_spai_v("MGBench mcsse | q1t | sort 0 | L4 | spai | V", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_spai_v("MGBench mcsse | q1t | sort 0 | L5 | spai | V", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_spai_v("MGBench mcsse | q1t | sort 0 | L6 | spai | V", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_spai_v("MGBench mcsse | q1t | sort 1 | L4 | spai | V", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_spai_v("MGBench mcsse | q1t | sort 1 | L5 | spai | V", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_spai_v("MGBench mcsse | q1t | sort 1 | L6 | spai | V", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_spai_v("MGBench mcsse | q1t | sort 2 | L4 | spai | V", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_spai_v("MGBench mcsse | q1t | sort 2 | L5 | spai | V", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_spai_v("MGBench mcsse | q1t | sort 2 | L6 | spai | V", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_spai_v("MGBench mcsse | q1t | sort 3 | L4 | spai | V", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_spai_v("MGBench mcsse | q1t | sort 3 | L5 | spai | V", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_spai_v("MGBench mcsse | q1t | sort 3 | L6 | spai | V", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_spai_v("MGBench mcsse | q1t | sort 4 | L4 | spai | V", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_spai_v("MGBench mcsse | q1t | sort 4 | L5 | spai | V", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_spai_v("MGBench mcsse | q1t | sort 4 | L6 | spai | V", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_sainv_v("MGBench mcsse | q1 | sort 0 | L5 | sainv | V", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_sainv_v("MGBench mcsse | q1 | sort 0 | L6 | sainv | V", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_v("MGBench mcsse | q1 | sort 0 | L7 | sainv | V", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_sainv_v("MGBench mcsse | q1 | sort 1 | L5 | sainv | V", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_sainv_v("MGBench mcsse | q1 | sort 1 | L6 | sainv | V", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_v("MGBench mcsse | q1 | sort 1 | L7 | sainv | V", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_sainv_v("MGBench mcsse | q1 | sort 2 | L5 | sainv | V", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_sainv_v("MGBench mcsse | q1 | sort 2 | L6 | sainv | V", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_sainv_v("MGBench mcsse | q1 | sort 2 | L7 | sainv | V", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_sainv_v("MGBench mcsse | q1 | sort 3 | L5 | sainv | V", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_sainv_v("MGBench mcsse | q1 | sort 3 | L6 | sainv | V", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_v("MGBench mcsse | q1 | sort 3 | L7 | sainv | V", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_sainv_v("MGBench mcsse | q1 | sort 4 | L5 | sainv | V", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_sainv_v("MGBench mcsse | q1 | sort 4 | L6 | sainv | V", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_sainv_v("MGBench mcsse | q1 | sort 4 | L7 | sainv | V", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_sainv_v("MGBench mcsse | q2 | sort 0 | L4 | sainv | V", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_sainv_v("MGBench mcsse | q2 | sort 0 | L5 | sainv | V", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_sainv_v("MGBench mcsse | q2 | sort 0 | L6 | sainv | V", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_sainv_v("MGBench mcsse | q2 | sort 1 | L4 | sainv | V", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_sainv_v("MGBench mcsse | q2 | sort 1 | L5 | sainv | V", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_sainv_v("MGBench mcsse | q2 | sort 1 | L6 | sainv | V", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_sainv_v("MGBench mcsse | q2 | sort 2 | L4 | sainv | V", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_sainv_v("MGBench mcsse | q2 | sort 2 | L5 | sainv | V", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_sainv_v("MGBench mcsse | q2 | sort 2 | L6 | sainv | V", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_sainv_v("MGBench mcsse | q2 | sort 3 | L4 | sainv | V", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_sainv_v("MGBench mcsse | q2 | sort 3 | L5 | sainv | V", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_sainv_v("MGBench mcsse | q2 | sort 3 | L6 | sainv | V", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_sainv_v("MGBench mcsse | q2 | sort 4 | L4 | sainv | V", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_sainv_v("MGBench mcsse | q2 | sort 4 | L5 | sainv | V", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_sainv_v("MGBench mcsse | q2 | sort 4 | L6 | sainv | V", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_sainv_v("MGBench mcsse | q1t | sort 0 | L4 | sainv | V", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_sainv_v("MGBench mcsse | q1t | sort 0 | L5 | sainv | V", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_sainv_v("MGBench mcsse | q1t | sort 0 | L6 | sainv | V", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_sainv_v("MGBench mcsse | q1t | sort 1 | L4 | sainv | V", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_sainv_v("MGBench mcsse | q1t | sort 1 | L5 | sainv | V", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_sainv_v("MGBench mcsse | q1t | sort 1 | L6 | sainv | V", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_sainv_v("MGBench mcsse | q1t | sort 2 | L4 | sainv | V", 3 , 2, 4, "sainv", 0.65);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_sainv_v("MGBench mcsse | q1t | sort 2 | L5 | sainv | V", 3 , 2, 5, "sainv", 0.65);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_sainv_v("MGBench mcsse | q1t | sort 2 | L6 | sainv | V", 3 , 2, 6, "sainv", 0.65);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_sainv_v("MGBench mcsse | q1t | sort 3 | L4 | sainv | V", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_sainv_v("MGBench mcsse | q1t | sort 3 | L5 | sainv | V", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_sainv_v("MGBench mcsse | q1t | sort 3 | L6 | sainv | V", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_sainv_v("MGBench mcsse | q1t | sort 4 | L4 | sainv | V", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_sainv_v("MGBench mcsse | q1t | sort 4 | L5 | sainv | V", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_sainv_v("MGBench mcsse | q1t | sort 4 | L6 | sainv | V", 3 , 4, 6, "sainv", 0.7);
#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l5_jac_v("MGBench cuda | q1 | sort 0 | L5 | jac | V", 1 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l6_jac_v("MGBench cuda | q1 | sort 0 | L6 | jac | V", 1 , 0, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_v("MGBench cuda | q1 | sort 0 | L7 | jac | V", 1 , 0, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l5_jac_v("MGBench cuda | q1 | sort 1 | L5 | jac | V", 1 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l6_jac_v("MGBench cuda | q1 | sort 1 | L6 | jac | V", 1 , 1, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_v("MGBench cuda | q1 | sort 1 | L7 | jac | V", 1 , 1, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort2_l5_jac_v("MGBench cuda | q1 | sort 2 | L5 | jac | V", 1 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort2_l6_jac_v("MGBench cuda | q1 | sort 2 | L6 | jac | V", 1 , 2, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort2_l7_jac_v("MGBench cuda | q1 | sort 2 | L7 | jac | V", 1 , 2, 7, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l5_jac_v("MGBench cuda | q1 | sort 3 | L5 | jac | V", 1 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l6_jac_v("MGBench cuda | q1 | sort 3 | L6 | jac | V", 1 , 3, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_v("MGBench cuda | q1 | sort 3 | L7 | jac | V", 1 , 3, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort4_l5_jac_v("MGBench cuda | q1 | sort 4 | L5 | jac | V", 1 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort4_l6_jac_v("MGBench cuda | q1 | sort 4 | L6 | jac | V", 1 , 4, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1_sort4_l7_jac_v("MGBench cuda | q1 | sort 4 | L7 | jac | V", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort0_l4_jac_v("MGBench cuda | q2 | sort 0 | L4 | jac | V", 2 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort0_l5_jac_v("MGBench cuda | q2 | sort 0 | L5 | jac | V", 2 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort0_l6_jac_v("MGBench cuda | q2 | sort 0 | L6 | jac | V", 2 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort1_l4_jac_v("MGBench cuda | q2 | sort 1 | L4 | jac | V", 2 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort1_l5_jac_v("MGBench cuda | q2 | sort 1 | L5 | jac | V", 2 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort1_l6_jac_v("MGBench cuda | q2 | sort 1 | L6 | jac | V", 2 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort2_l4_jac_v("MGBench cuda | q2 | sort 2 | L4 | jac | V", 2 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort2_l5_jac_v("MGBench cuda | q2 | sort 2 | L5 | jac | V", 2 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort2_l6_jac_v("MGBench cuda | q2 | sort 2 | L6 | jac | V", 2 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort3_l4_jac_v("MGBench cuda | q2 | sort 3 | L4 | jac | V", 2 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort3_l5_jac_v("MGBench cuda | q2 | sort 3 | L5 | jac | V", 2 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort3_l6_jac_v("MGBench cuda | q2 | sort 3 | L6 | jac | V", 2 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort4_l4_jac_v("MGBench cuda | q2 | sort 4 | L4 | jac | V", 2 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort4_l5_jac_v("MGBench cuda | q2 | sort 4 | L5 | jac | V", 2 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q2_sort4_l6_jac_v("MGBench cuda | q2 | sort 4 | L6 | jac | V", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l4_jac_v("MGBench cuda | q1t | sort 0 | L4 | jac | V", 3 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l5_jac_v("MGBench cuda | q1t | sort 0 | L5 | jac | V", 3 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l6_jac_v("MGBench cuda | q1t | sort 0 | L6 | jac | V", 3 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l4_jac_v("MGBench cuda | q1t | sort 1 | L4 | jac | V", 3 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l5_jac_v("MGBench cuda | q1t | sort 1 | L5 | jac | V", 3 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l6_jac_v("MGBench cuda | q1t | sort 1 | L6 | jac | V", 3 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l4_jac_v("MGBench cuda | q1t | sort 2 | L4 | jac | V", 3 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l5_jac_v("MGBench cuda | q1t | sort 2 | L5 | jac | V", 3 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l6_jac_v("MGBench cuda | q1t | sort 2 | L6 | jac | V", 3 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l4_jac_v("MGBench cuda | q1t | sort 3 | L4 | jac | V", 3 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l5_jac_v("MGBench cuda | q1t | sort 3 | L5 | jac | V", 3 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l6_jac_v("MGBench cuda | q1t | sort 3 | L6 | jac | V", 3 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l4_jac_v("MGBench cuda | q1t | sort 4 | L4 | jac | V", 3 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l5_jac_v("MGBench cuda | q1t | sort 4 | L5 | jac | V", 3 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l6_jac_v("MGBench cuda | q1t | sort 4 | L6 | jac | V", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_spai_v("MGBench cuda | q1 | sort 0 | L5 | spai | V", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_spai_v("MGBench cuda | q1 | sort 0 | L6 | spai | V", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_v("MGBench cuda | q1 | sort 0 | L7 | spai | V", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_spai_v("MGBench cuda | q1 | sort 1 | L5 | spai | V", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_spai_v("MGBench cuda | q1 | sort 1 | L6 | spai | V", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_v("MGBench cuda | q1 | sort 1 | L7 | spai | V", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_spai_v("MGBench cuda | q1 | sort 2 | L5 | spai | V", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_spai_v("MGBench cuda | q1 | sort 2 | L6 | spai | V", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_spai_v("MGBench cuda | q1 | sort 2 | L7 | spai | V", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_spai_v("MGBench cuda | q1 | sort 3 | L5 | spai | V", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_spai_v("MGBench cuda | q1 | sort 3 | L6 | spai | V", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_v("MGBench cuda | q1 | sort 3 | L7 | spai | V", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_spai_v("MGBench cuda | q1 | sort 4 | L5 | spai | V", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_spai_v("MGBench cuda | q1 | sort 4 | L6 | spai | V", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_spai_v("MGBench cuda | q1 | sort 4 | L7 | spai | V", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_spai_v("MGBench cuda | q2 | sort 0 | L4 | spai | V", 2 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_spai_v("MGBench cuda | q2 | sort 0 | L5 | spai | V", 2 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_spai_v("MGBench cuda | q2 | sort 0 | L6 | spai | V", 2 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_spai_v("MGBench cuda | q2 | sort 1 | L4 | spai | V", 2 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_spai_v("MGBench cuda | q2 | sort 1 | L5 | spai | V", 2 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_spai_v("MGBench cuda | q2 | sort 1 | L6 | spai | V", 2 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_spai_v("MGBench cuda | q2 | sort 2 | L4 | spai | V", 2 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_spai_v("MGBench cuda | q2 | sort 2 | L5 | spai | V", 2 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_spai_v("MGBench cuda | q2 | sort 2 | L6 | spai | V", 2 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_spai_v("MGBench cuda | q2 | sort 3 | L4 | spai | V", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_spai_v("MGBench cuda | q2 | sort 3 | L5 | spai | V", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_spai_v("MGBench cuda | q2 | sort 3 | L6 | spai | V", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_spai_v("MGBench cuda | q2 | sort 4 | L4 | spai | V", 2 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_spai_v("MGBench cuda | q2 | sort 4 | L5 | spai | V", 2 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_spai_v("MGBench cuda | q2 | sort 4 | L6 | spai | V", 2 , 4, 6, "spai_grote", 1.);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_spai_v("MGBench cuda | q1t | sort 0 | L4 | spai | V", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_spai_v("MGBench cuda | q1t | sort 0 | L5 | spai | V", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_spai_v("MGBench cuda | q1t | sort 0 | L6 | spai | V", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_spai_v("MGBench cuda | q1t | sort 1 | L4 | spai | V", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_spai_v("MGBench cuda | q1t | sort 1 | L5 | spai | V", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_spai_v("MGBench cuda | q1t | sort 1 | L6 | spai | V", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_spai_v("MGBench cuda | q1t | sort 2 | L4 | spai | V", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_spai_v("MGBench cuda | q1t | sort 2 | L5 | spai | V", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_spai_v("MGBench cuda | q1t | sort 2 | L6 | spai | V", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_spai_v("MGBench cuda | q1t | sort 3 | L4 | spai | V", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_spai_v("MGBench cuda | q1t | sort 3 | L5 | spai | V", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_spai_v("MGBench cuda | q1t | sort 3 | L6 | spai | V", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_spai_v("MGBench cuda | q1t | sort 4 | L4 | spai | V", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_spai_v("MGBench cuda | q1t | sort 4 | L5 | spai | V", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_spai_v("MGBench cuda | q1t | sort 4 | L6 | spai | V", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_sainv_v("MGBench cuda | q1 | sort 0 | L5 | sainv | V", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_sainv_v("MGBench cuda | q1 | sort 0 | L6 | sainv | V", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_v("MGBench cuda | q1 | sort 0 | L7 | sainv | V", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_sainv_v("MGBench cuda | q1 | sort 1 | L5 | sainv | V", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_sainv_v("MGBench cuda | q1 | sort 1 | L6 | sainv | V", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_v("MGBench cuda | q1 | sort 1 | L7 | sainv | V", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_sainv_v("MGBench cuda | q1 | sort 2 | L5 | sainv | V", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_sainv_v("MGBench cuda | q1 | sort 2 | L6 | sainv | V", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_sainv_v("MGBench cuda | q1 | sort 2 | L7 | sainv | V", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_sainv_v("MGBench cuda | q1 | sort 3 | L5 | sainv | V", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_sainv_v("MGBench cuda | q1 | sort 3 | L6 | sainv | V", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_v("MGBench cuda | q1 | sort 3 | L7 | sainv | V", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_sainv_v("MGBench cuda | q1 | sort 4 | L5 | sainv | V", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_sainv_v("MGBench cuda | q1 | sort 4 | L6 | sainv | V", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_sainv_v("MGBench cuda | q1 | sort 4 | L7 | sainv | V", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_sainv_v("MGBench cuda | q2 | sort 0 | L4 | sainv | V", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_sainv_v("MGBench cuda | q2 | sort 0 | L5 | sainv | V", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_sainv_v("MGBench cuda | q2 | sort 0 | L6 | sainv | V", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_sainv_v("MGBench cuda | q2 | sort 1 | L4 | sainv | V", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_sainv_v("MGBench cuda | q2 | sort 1 | L5 | sainv | V", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_sainv_v("MGBench cuda | q2 | sort 1 | L6 | sainv | V", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_sainv_v("MGBench cuda | q2 | sort 2 | L4 | sainv | V", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_sainv_v("MGBench cuda | q2 | sort 2 | L5 | sainv | V", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_sainv_v("MGBench cuda | q2 | sort 2 | L6 | sainv | V", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_sainv_v("MGBench cuda | q2 | sort 3 | L4 | sainv | V", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_sainv_v("MGBench cuda | q2 | sort 3 | L5 | sainv | V", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_sainv_v("MGBench cuda | q2 | sort 3 | L6 | sainv | V", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_sainv_v("MGBench cuda | q2 | sort 4 | L4 | sainv | V", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_sainv_v("MGBench cuda | q2 | sort 4 | L5 | sainv | V", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_sainv_v("MGBench cuda | q2 | sort 4 | L6 | sainv | V", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_sainv_v("MGBench cuda | q1t | sort 0 | L4 | sainv | V", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_sainv_v("MGBench cuda | q1t | sort 0 | L5 | sainv | V", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_sainv_v("MGBench cuda | q1t | sort 0 | L6 | sainv | V", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_sainv_v("MGBench cuda | q1t | sort 1 | L4 | sainv | V", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_sainv_v("MGBench cuda | q1t | sort 1 | L5 | sainv | V", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_sainv_v("MGBench cuda | q1t | sort 1 | L6 | sainv | V", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_sainv_v("MGBench cuda | q1t | sort 2 | L4 | sainv | V", 3 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_sainv_v("MGBench cuda | q1t | sort 2 | L5 | sainv | V", 3 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_sainv_v("MGBench cuda | q1t | sort 2 | L6 | sainv | V", 3 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_sainv_v("MGBench cuda | q1t | sort 3 | L4 | sainv | V", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_sainv_v("MGBench cuda | q1t | sort 3 | L5 | sainv | V", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_sainv_v("MGBench cuda | q1t | sort 3 | L6 | sainv | V", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_sainv_v("MGBench cuda | q1t | sort 4 | L4 | sainv | V", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_sainv_v("MGBench cuda | q1t | sort 4 | L5 | sainv | V", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_sainv_v("MGBench cuda | q1t | sort 4 | L6 | sainv | V", 3 , 4, 6, "sainv", 0.7);
#endif

#ifdef HONEI_SSE
//JAC
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l5_jac_w("MGBench mcsse | q1 | sort 0 | L5 | jac | W", 1 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l6_jac_w("MGBench mcsse | q1 | sort 0 | L6 | jac | W", 1 , 0, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_w("MGBench mcsse | q1 | sort 0 | L7 | jac | W", 1 , 0, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l5_jac_w("MGBench mcsse | q1 | sort 1 | L5 | jac | W", 1 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l6_jac_w("MGBench mcsse | q1 | sort 1 | L6 | jac | W", 1 , 1, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_w("MGBench mcsse | q1 | sort 1 | L7 | jac | W", 1 , 1, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort2_l5_jac_w("MGBench mcsse | q1 | sort 2 | L5 | jac | W", 1 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort2_l6_jac_w("MGBench mcsse | q1 | sort 2 | L6 | jac | W", 1 , 2, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort2_l7_jac_w("MGBench mcsse | q1 | sort 2 | L7 | jac | W", 1 , 2, 7, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l5_jac_w("MGBench mcsse | q1 | sort 3 | L5 | jac | W", 1 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l6_jac_w("MGBench mcsse | q1 | sort 3 | L6 | jac | W", 1 , 3, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_w("MGBench mcsse | q1 | sort 3 | L7 | jac | W", 1 , 3, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort4_l5_jac_w("MGBench mcsse | q1 | sort 4 | L5 | jac | W", 1 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort4_l6_jac_w("MGBench mcsse | q1 | sort 4 | L6 | jac | W", 1 , 4, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1_sort4_l7_jac_w("MGBench mcsse | q1 | sort 4 | L7 | jac | W", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort0_l4_jac_w("MGBench mcsse | q2 | sort 0 | L4 | jac | W", 2 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort0_l5_jac_w("MGBench mcsse | q2 | sort 0 | L5 | jac | W", 2 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort0_l6_jac_w("MGBench mcsse | q2 | sort 0 | L6 | jac | W", 2 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort1_l4_jac_w("MGBench mcsse | q2 | sort 1 | L4 | jac | W", 2 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort1_l5_jac_w("MGBench mcsse | q2 | sort 1 | L5 | jac | W", 2 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort1_l6_jac_w("MGBench mcsse | q2 | sort 1 | L6 | jac | W", 2 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort2_l4_jac_w("MGBench mcsse | q2 | sort 2 | L4 | jac | W", 2 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort2_l5_jac_w("MGBench mcsse | q2 | sort 2 | L5 | jac | W", 2 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort2_l6_jac_w("MGBench mcsse | q2 | sort 2 | L6 | jac | W", 2 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort3_l4_jac_w("MGBench mcsse | q2 | sort 3 | L4 | jac | W", 2 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort3_l5_jac_w("MGBench mcsse | q2 | sort 3 | L5 | jac | W", 2 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort3_l6_jac_w("MGBench mcsse | q2 | sort 3 | L6 | jac | W", 2 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort4_l4_jac_w("MGBench mcsse | q2 | sort 4 | L4 | jac | W", 2 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort4_l5_jac_w("MGBench mcsse | q2 | sort 4 | L5 | jac | W", 2 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q2_sort4_l6_jac_w("MGBench mcsse | q2 | sort 4 | L6 | jac | W", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort0_l4_jac_w("MGBench mcsse | q1t | sort 0 | L4 | jac | W", 3 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort0_l5_jac_w("MGBench mcsse | q1t | sort 0 | L5 | jac | W", 3 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort0_l6_jac_w("MGBench mcsse | q1t | sort 0 | L6 | jac | W", 3 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort1_l4_jac_w("MGBench mcsse | q1t | sort 1 | L4 | jac | W", 3 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort1_l5_jac_w("MGBench mcsse | q1t | sort 1 | L5 | jac | W", 3 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort1_l6_jac_w("MGBench mcsse | q1t | sort 1 | L6 | jac | W", 3 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort2_l4_jac_w("MGBench mcsse | q1t | sort 2 | L4 | jac | W", 3 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort2_l5_jac_w("MGBench mcsse | q1t | sort 2 | L5 | jac | W", 3 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort2_l6_jac_w("MGBench mcsse | q1t | sort 2 | L6 | jac | W", 3 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort3_l4_jac_w("MGBench mcsse | q1t | sort 3 | L4 | jac | W", 3 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort3_l5_jac_w("MGBench mcsse | q1t | sort 3 | L5 | jac | W", 3 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort3_l6_jac_w("MGBench mcsse | q1t | sort 3 | L6 | jac | W", 3 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort4_l4_jac_w("MGBench mcsse | q1t | sort 4 | L4 | jac | W", 3 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort4_l5_jac_w("MGBench mcsse | q1t | sort 4 | L5 | jac | W", 3 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, DenseVector<double> > mcsse_q1t_sort4_l6_jac_w("MGBench mcsse | q1t | sort 4 | L6 | jac | W", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_spai_w("MGBench mcsse | q1 | sort 0 | L5 | spai | W", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_spai_w("MGBench mcsse | q1 | sort 0 | L6 | spai | W", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_w("MGBench mcsse | q1 | sort 0 | L7 | spai | W", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_spai_w("MGBench mcsse | q1 | sort 1 | L5 | spai | W", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_spai_w("MGBench mcsse | q1 | sort 1 | L6 | spai | W", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_w("MGBench mcsse | q1 | sort 1 | L7 | spai | W", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_spai_w("MGBench mcsse | q1 | sort 2 | L5 | spai | W", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_spai_w("MGBench mcsse | q1 | sort 2 | L6 | spai | W", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_spai_w("MGBench mcsse | q1 | sort 2 | L7 | spai | W", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_spai_w("MGBench mcsse | q1 | sort 3 | L5 | spai | W", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_spai_w("MGBench mcsse | q1 | sort 3 | L6 | spai | W", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_w("MGBench mcsse | q1 | sort 3 | L7 | spai | W", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_spai_w("MGBench mcsse | q1 | sort 4 | L5 | spai | W", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_spai_w("MGBench mcsse | q1 | sort 4 | L6 | spai | W", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_spai_w("MGBench mcsse | q1 | sort 4 | L7 | spai | W", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_spai_w("MGBench mcsse | q2 | sort 0 | L4 | spai | W", 2 , 0, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_spai_w("MGBench mcsse | q2 | sort 0 | L5 | spai | W", 2 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_spai_w("MGBench mcsse | q2 | sort 0 | L6 | spai | W", 2 , 0, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_spai_w("MGBench mcsse | q2 | sort 1 | L4 | spai | W", 2 , 1, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_spai_w("MGBench mcsse | q2 | sort 1 | L5 | spai | W", 2 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_spai_w("MGBench mcsse | q2 | sort 1 | L6 | spai | W", 2 , 1, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_spai_w("MGBench mcsse | q2 | sort 2 | L4 | spai | W", 2 , 2, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_spai_w("MGBench mcsse | q2 | sort 2 | L5 | spai | W", 2 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_spai_w("MGBench mcsse | q2 | sort 2 | L6 | spai | W", 2 , 2, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_spai_w("MGBench mcsse | q2 | sort 3 | L4 | spai | W", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_spai_w("MGBench mcsse | q2 | sort 3 | L5 | spai | W", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_spai_w("MGBench mcsse | q2 | sort 3 | L6 | spai | W", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_spai_w("MGBench mcsse | q2 | sort 4 | L4 | spai | W", 2 , 4, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_spai_w("MGBench mcsse | q2 | sort 4 | L5 | spai | W", 2 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_spai_w("MGBench mcsse | q2 | sort 4 | L6 | spai | W", 2 , 4, 6, "spai_grote", 1.);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_spai_w("MGBench mcsse | q1t | sort 0 | L4 | spai | W", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_spai_w("MGBench mcsse | q1t | sort 0 | L5 | spai | W", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_spai_w("MGBench mcsse | q1t | sort 0 | L6 | spai | W", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_spai_w("MGBench mcsse | q1t | sort 1 | L4 | spai | W", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_spai_w("MGBench mcsse | q1t | sort 1 | L5 | spai | W", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_spai_w("MGBench mcsse | q1t | sort 1 | L6 | spai | W", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_spai_w("MGBench mcsse | q1t | sort 2 | L4 | spai | W", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_spai_w("MGBench mcsse | q1t | sort 2 | L5 | spai | W", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_spai_w("MGBench mcsse | q1t | sort 2 | L6 | spai | W", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_spai_w("MGBench mcsse | q1t | sort 3 | L4 | spai | W", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_spai_w("MGBench mcsse | q1t | sort 3 | L5 | spai | W", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_spai_w("MGBench mcsse | q1t | sort 3 | L6 | spai | W", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_spai_w("MGBench mcsse | q1t | sort 4 | L4 | spai | W", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_spai_w("MGBench mcsse | q1t | sort 4 | L5 | spai | W", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_spai_w("MGBench mcsse | q1t | sort 4 | L6 | spai | W", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_sainv_w("MGBench mcsse | q1 | sort 0 | L5 | sainv | W", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_sainv_w("MGBench mcsse | q1 | sort 0 | L6 | sainv | W", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_w("MGBench mcsse | q1 | sort 0 | L7 | sainv | W", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_sainv_w("MGBench mcsse | q1 | sort 1 | L5 | sainv | W", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_sainv_w("MGBench mcsse | q1 | sort 1 | L6 | sainv | W", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_w("MGBench mcsse | q1 | sort 1 | L7 | sainv | W", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_sainv_w("MGBench mcsse | q1 | sort 2 | L5 | sainv | W", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_sainv_w("MGBench mcsse | q1 | sort 2 | L6 | sainv | W", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_sainv_w("MGBench mcsse | q1 | sort 2 | L7 | sainv | W", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_sainv_w("MGBench mcsse | q1 | sort 3 | L5 | sainv | W", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_sainv_w("MGBench mcsse | q1 | sort 3 | L6 | sainv | W", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_w("MGBench mcsse | q1 | sort 3 | L7 | sainv | W", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_sainv_w("MGBench mcsse | q1 | sort 4 | L5 | sainv | W", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_sainv_w("MGBench mcsse | q1 | sort 4 | L6 | sainv | W", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_sainv_w("MGBench mcsse | q1 | sort 4 | L7 | sainv | W", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_sainv_w("MGBench mcsse | q2 | sort 0 | L4 | sainv | W", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_sainv_w("MGBench mcsse | q2 | sort 0 | L5 | sainv | W", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_sainv_w("MGBench mcsse | q2 | sort 0 | L6 | sainv | W", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_sainv_w("MGBench mcsse | q2 | sort 1 | L4 | sainv | W", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_sainv_w("MGBench mcsse | q2 | sort 1 | L5 | sainv | W", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_sainv_w("MGBench mcsse | q2 | sort 1 | L6 | sainv | W", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_sainv_w("MGBench mcsse | q2 | sort 2 | L4 | sainv | W", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_sainv_w("MGBench mcsse | q2 | sort 2 | L5 | sainv | W", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_sainv_w("MGBench mcsse | q2 | sort 2 | L6 | sainv | W", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_sainv_w("MGBench mcsse | q2 | sort 3 | L4 | sainv | W", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_sainv_w("MGBench mcsse | q2 | sort 3 | L5 | sainv | W", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_sainv_w("MGBench mcsse | q2 | sort 3 | L6 | sainv | W", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_sainv_w("MGBench mcsse | q2 | sort 4 | L4 | sainv | W", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_sainv_w("MGBench mcsse | q2 | sort 4 | L5 | sainv | W", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_sainv_w("MGBench mcsse | q2 | sort 4 | L6 | sainv | W", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_sainv_w("MGBench mcsse | q1t | sort 0 | L4 | sainv | W", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_sainv_w("MGBench mcsse | q1t | sort 0 | L5 | sainv | W", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_sainv_w("MGBench mcsse | q1t | sort 0 | L6 | sainv | W", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_sainv_w("MGBench mcsse | q1t | sort 1 | L4 | sainv | W", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_sainv_w("MGBench mcsse | q1t | sort 1 | L5 | sainv | W", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_sainv_w("MGBench mcsse | q1t | sort 1 | L6 | sainv | W", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_sainv_w("MGBench mcsse | q1t | sort 2 | L4 | sainv | W", 3 , 2, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_sainv_w("MGBench mcsse | q1t | sort 2 | L5 | sainv | W", 3 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_sainv_w("MGBench mcsse | q1t | sort 2 | L6 | sainv | W", 3 , 2, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_sainv_w("MGBench mcsse | q1t | sort 3 | L4 | sainv | W", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_sainv_w("MGBench mcsse | q1t | sort 3 | L5 | sainv | W", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_sainv_w("MGBench mcsse | q1t | sort 3 | L6 | sainv | W", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_sainv_w("MGBench mcsse | q1t | sort 4 | L4 | sainv | W", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_sainv_w("MGBench mcsse | q1t | sort 4 | L5 | sainv | W", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_sainv_w("MGBench mcsse | q1t | sort 4 | L6 | sainv | W", 3 , 4, 6, "sainv", 0.7);
#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l5_jac_w("MGBench cuda | q1 | sort 0 | L5 | jac | W", 1 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l6_jac_w("MGBench cuda | q1 | sort 0 | L6 | jac | W", 1 , 0, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_w("MGBench cuda | q1 | sort 0 | L7 | jac | W", 1 , 0, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l5_jac_w("MGBench cuda | q1 | sort 1 | L5 | jac | W", 1 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l6_jac_w("MGBench cuda | q1 | sort 1 | L6 | jac | W", 1 , 1, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_w("MGBench cuda | q1 | sort 1 | L7 | jac | W", 1 , 1, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort2_l5_jac_w("MGBench cuda | q1 | sort 2 | L5 | jac | W", 1 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort2_l6_jac_w("MGBench cuda | q1 | sort 2 | L6 | jac | W", 1 , 2, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort2_l7_jac_w("MGBench cuda | q1 | sort 2 | L7 | jac | W", 1 , 2, 7, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l5_jac_w("MGBench cuda | q1 | sort 3 | L5 | jac | W", 1 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l6_jac_w("MGBench cuda | q1 | sort 3 | L6 | jac | W", 1 , 3, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_w("MGBench cuda | q1 | sort 3 | L7 | jac | W", 1 , 3, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort4_l5_jac_w("MGBench cuda | q1 | sort 4 | L5 | jac | W", 1 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort4_l6_jac_w("MGBench cuda | q1 | sort 4 | L6 | jac | W", 1 , 4, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1_sort4_l7_jac_w("MGBench cuda | q1 | sort 4 | L7 | jac | W", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort0_l4_jac_w("MGBench cuda | q2 | sort 0 | L4 | jac | W", 2 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort0_l5_jac_w("MGBench cuda | q2 | sort 0 | L5 | jac | W", 2 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort0_l6_jac_w("MGBench cuda | q2 | sort 0 | L6 | jac | W", 2 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort1_l4_jac_w("MGBench cuda | q2 | sort 1 | L4 | jac | W", 2 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort1_l5_jac_w("MGBench cuda | q2 | sort 1 | L5 | jac | W", 2 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort1_l6_jac_w("MGBench cuda | q2 | sort 1 | L6 | jac | W", 2 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort2_l4_jac_w("MGBench cuda | q2 | sort 2 | L4 | jac | W", 2 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort2_l5_jac_w("MGBench cuda | q2 | sort 2 | L5 | jac | W", 2 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort2_l6_jac_w("MGBench cuda | q2 | sort 2 | L6 | jac | W", 2 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort3_l4_jac_w("MGBench cuda | q2 | sort 3 | L4 | jac | W", 2 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort3_l5_jac_w("MGBench cuda | q2 | sort 3 | L5 | jac | W", 2 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort3_l6_jac_w("MGBench cuda | q2 | sort 3 | L6 | jac | W", 2 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort4_l4_jac_w("MGBench cuda | q2 | sort 4 | L4 | jac | W", 2 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort4_l5_jac_w("MGBench cuda | q2 | sort 4 | L5 | jac | W", 2 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q2_sort4_l6_jac_w("MGBench cuda | q2 | sort 4 | L6 | jac | W", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort0_l4_jac_w("MGBench cuda | q1t | sort 0 | L4 | jac | W", 3 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort0_l5_jac_w("MGBench cuda | q1t | sort 0 | L5 | jac | W", 3 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort0_l6_jac_w("MGBench cuda | q1t | sort 0 | L6 | jac | W", 3 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort1_l4_jac_w("MGBench cuda | q1t | sort 1 | L4 | jac | W", 3 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort1_l5_jac_w("MGBench cuda | q1t | sort 1 | L5 | jac | W", 3 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort1_l6_jac_w("MGBench cuda | q1t | sort 1 | L6 | jac | W", 3 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort2_l4_jac_w("MGBench cuda | q1t | sort 2 | L4 | jac | W", 3 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort2_l5_jac_w("MGBench cuda | q1t | sort 2 | L5 | jac | W", 3 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort2_l6_jac_w("MGBench cuda | q1t | sort 2 | L6 | jac | W", 3 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort3_l4_jac_w("MGBench cuda | q1t | sort 3 | L4 | jac | W", 3 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort3_l5_jac_w("MGBench cuda | q1t | sort 3 | L5 | jac | W", 3 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort3_l6_jac_w("MGBench cuda | q1t | sort 3 | L6 | jac | W", 3 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort4_l4_jac_w("MGBench cuda | q1t | sort 4 | L4 | jac | W", 3 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort4_l5_jac_w("MGBench cuda | q1t | sort 4 | L5 | jac | W", 3 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, DenseVector<double> > cuda_q1t_sort4_l6_jac_w("MGBench cuda | q1t | sort 4 | L6 | jac | W", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_spai_w("MGBench cuda | q1 | sort 0 | L5 | spai | W", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_spai_w("MGBench cuda | q1 | sort 0 | L6 | spai | W", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_w("MGBench cuda | q1 | sort 0 | L7 | spai | W", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_spai_w("MGBench cuda | q1 | sort 1 | L5 | spai | W", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_spai_w("MGBench cuda | q1 | sort 1 | L6 | spai | W", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_w("MGBench cuda | q1 | sort 1 | L7 | spai | W", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_spai_w("MGBench cuda | q1 | sort 2 | L5 | spai | W", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_spai_w("MGBench cuda | q1 | sort 2 | L6 | spai | W", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_spai_w("MGBench cuda | q1 | sort 2 | L7 | spai | W", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_spai_w("MGBench cuda | q1 | sort 3 | L5 | spai | W", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_spai_w("MGBench cuda | q1 | sort 3 | L6 | spai | W", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_w("MGBench cuda | q1 | sort 3 | L7 | spai | W", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_spai_w("MGBench cuda | q1 | sort 4 | L5 | spai | W", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_spai_w("MGBench cuda | q1 | sort 4 | L6 | spai | W", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_spai_w("MGBench cuda | q1 | sort 4 | L7 | spai | W", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_spai_w("MGBench cuda | q2 | sort 0 | L4 | spai | W", 2 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_spai_w("MGBench cuda | q2 | sort 0 | L5 | spai | W", 2 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_spai_w("MGBench cuda | q2 | sort 0 | L6 | spai | W", 2 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_spai_w("MGBench cuda | q2 | sort 1 | L4 | spai | W", 2 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_spai_w("MGBench cuda | q2 | sort 1 | L5 | spai | W", 2 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_spai_w("MGBench cuda | q2 | sort 1 | L6 | spai | W", 2 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_spai_w("MGBench cuda | q2 | sort 2 | L4 | spai | W", 2 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_spai_w("MGBench cuda | q2 | sort 2 | L5 | spai | W", 2 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_spai_w("MGBench cuda | q2 | sort 2 | L6 | spai | W", 2 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_spai_w("MGBench cuda | q2 | sort 3 | L4 | spai | W", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_spai_w("MGBench cuda | q2 | sort 3 | L5 | spai | W", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_spai_w("MGBench cuda | q2 | sort 3 | L6 | spai | W", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_spai_w("MGBench cuda | q2 | sort 4 | L4 | spai | W", 2 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_spai_w("MGBench cuda | q2 | sort 4 | L5 | spai | W", 2 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_spai_w("MGBench cuda | q2 | sort 4 | L6 | spai | W", 2 , 4, 6, "spai_grote", 1.);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_spai_w("MGBench cuda | q1t | sort 0 | L4 | spai | W", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_spai_w("MGBench cuda | q1t | sort 0 | L5 | spai | W", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_spai_w("MGBench cuda | q1t | sort 0 | L6 | spai | W", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_spai_w("MGBench cuda | q1t | sort 1 | L4 | spai | W", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_spai_w("MGBench cuda | q1t | sort 1 | L5 | spai | W", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_spai_w("MGBench cuda | q1t | sort 1 | L6 | spai | W", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_spai_w("MGBench cuda | q1t | sort 2 | L4 | spai | W", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_spai_w("MGBench cuda | q1t | sort 2 | L5 | spai | W", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_spai_w("MGBench cuda | q1t | sort 2 | L6 | spai | W", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_spai_w("MGBench cuda | q1t | sort 3 | L4 | spai | W", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_spai_w("MGBench cuda | q1t | sort 3 | L5 | spai | W", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_spai_w("MGBench cuda | q1t | sort 3 | L6 | spai | W", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_spai_w("MGBench cuda | q1t | sort 4 | L4 | spai | W", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_spai_w("MGBench cuda | q1t | sort 4 | L5 | spai | W", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_spai_w("MGBench cuda | q1t | sort 4 | L6 | spai | W", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_sainv_w("MGBench cuda | q1 | sort 0 | L5 | sainv | W", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_sainv_w("MGBench cuda | q1 | sort 0 | L6 | sainv | W", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_w("MGBench cuda | q1 | sort 0 | L7 | sainv | W", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_sainv_w("MGBench cuda | q1 | sort 1 | L5 | sainv | W", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_sainv_w("MGBench cuda | q1 | sort 1 | L6 | sainv | W", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_w("MGBench cuda | q1 | sort 1 | L7 | sainv | W", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_sainv_w("MGBench cuda | q1 | sort 2 | L5 | sainv | W", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_sainv_w("MGBench cuda | q1 | sort 2 | L6 | sainv | W", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_sainv_w("MGBench cuda | q1 | sort 2 | L7 | sainv | W", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_sainv_w("MGBench cuda | q1 | sort 3 | L5 | sainv | W", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_sainv_w("MGBench cuda | q1 | sort 3 | L6 | sainv | W", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_w("MGBench cuda | q1 | sort 3 | L7 | sainv | W", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_sainv_w("MGBench cuda | q1 | sort 4 | L5 | sainv | W", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_sainv_w("MGBench cuda | q1 | sort 4 | L6 | sainv | W", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_sainv_w("MGBench cuda | q1 | sort 4 | L7 | sainv | W", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_sainv_w("MGBench cuda | q2 | sort 0 | L4 | sainv | W", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_sainv_w("MGBench cuda | q2 | sort 0 | L5 | sainv | W", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_sainv_w("MGBench cuda | q2 | sort 0 | L6 | sainv | W", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_sainv_w("MGBench cuda | q2 | sort 1 | L4 | sainv | W", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_sainv_w("MGBench cuda | q2 | sort 1 | L5 | sainv | W", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_sainv_w("MGBench cuda | q2 | sort 1 | L6 | sainv | W", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_sainv_w("MGBench cuda | q2 | sort 2 | L4 | sainv | W", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_sainv_w("MGBench cuda | q2 | sort 2 | L5 | sainv | W", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_sainv_w("MGBench cuda | q2 | sort 2 | L6 | sainv | W", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_sainv_w("MGBench cuda | q2 | sort 3 | L4 | sainv | W", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_sainv_w("MGBench cuda | q2 | sort 3 | L5 | sainv | W", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_sainv_w("MGBench cuda | q2 | sort 3 | L6 | sainv | W", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_sainv_w("MGBench cuda | q2 | sort 4 | L4 | sainv | W", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_sainv_w("MGBench cuda | q2 | sort 4 | L5 | sainv | W", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_sainv_w("MGBench cuda | q2 | sort 4 | L6 | sainv | W", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_sainv_w("MGBench cuda | q1t | sort 0 | L4 | sainv | W", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_sainv_w("MGBench cuda | q1t | sort 0 | L5 | sainv | W", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_sainv_w("MGBench cuda | q1t | sort 0 | L6 | sainv | W", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_sainv_w("MGBench cuda | q1t | sort 1 | L4 | sainv | W", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_sainv_w("MGBench cuda | q1t | sort 1 | L5 | sainv | W", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_sainv_w("MGBench cuda | q1t | sort 1 | L6 | sainv | W", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_sainv_w("MGBench cuda | q1t | sort 2 | L4 | sainv | W", 3 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_sainv_w("MGBench cuda | q1t | sort 2 | L5 | sainv | W", 3 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_sainv_w("MGBench cuda | q1t | sort 2 | L6 | sainv | W", 3 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_sainv_w("MGBench cuda | q1t | sort 3 | L4 | sainv | W", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_sainv_w("MGBench cuda | q1t | sort 3 | L5 | sainv | W", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_sainv_w("MGBench cuda | q1t | sort 3 | L6 | sainv | W", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_sainv_w("MGBench cuda | q1t | sort 4 | L4 | sainv | W", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_sainv_w("MGBench cuda | q1t | sort 4 | L5 | sainv | W", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::W::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_sainv_w("MGBench cuda | q1t | sort 4 | L6 | sainv | W", 3 , 4, 6, "sainv", 0.7);
#endif

#ifdef HONEI_SSE
//JAC
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l5_jac_f("MGBench mcsse | q1 | sort 0 | L5 | jac | F", 1 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l6_jac_f("MGBench mcsse | q1 | sort 0 | L6 | jac | F", 1 , 0, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort0_l7_jac_f("MGBench mcsse | q1 | sort 0 | L7 | jac | F", 1 , 0, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l5_jac_f("MGBench mcsse | q1 | sort 1 | L5 | jac | F", 1 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l6_jac_f("MGBench mcsse | q1 | sort 1 | L6 | jac | F", 1 , 1, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort1_l7_jac_f("MGBench mcsse | q1 | sort 1 | L7 | jac | F", 1 , 1, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l5_jac_f("MGBench mcsse | q1 | sort 2 | L5 | jac | F", 1 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l6_jac_f("MGBench mcsse | q1 | sort 2 | L6 | jac | F", 1 , 2, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l7_jac_f("MGBench mcsse | q1 | sort 2 | L7 | jac | F", 1 , 2, 7, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l5_jac_f("MGBench mcsse | q1 | sort 3 | L5 | jac | F", 1 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l6_jac_f("MGBench mcsse | q1 | sort 3 | L6 | jac | F", 1 , 3, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort3_l7_jac_f("MGBench mcsse | q1 | sort 3 | L7 | jac | F", 1 , 3, 7, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l5_jac_f("MGBench mcsse | q1 | sort 4 | L5 | jac | F", 1 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l6_jac_f("MGBench mcsse | q1 | sort 4 | L6 | jac | F", 1 , 4, 6, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1_sort4_l7_jac_f("MGBench mcsse | q1 | sort 4 | L7 | jac | F", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l4_jac_f("MGBench mcsse | q2 | sort 0 | L4 | jac | F", 2 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l5_jac_f("MGBench mcsse | q2 | sort 0 | L5 | jac | F", 2 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort0_l6_jac_f("MGBench mcsse | q2 | sort 0 | L6 | jac | F", 2 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l4_jac_f("MGBench mcsse | q2 | sort 1 | L4 | jac | F", 2 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l5_jac_f("MGBench mcsse | q2 | sort 1 | L5 | jac | F", 2 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort1_l6_jac_f("MGBench mcsse | q2 | sort 1 | L6 | jac | F", 2 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l4_jac_f("MGBench mcsse | q2 | sort 2 | L4 | jac | F", 2 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l5_jac_f("MGBench mcsse | q2 | sort 2 | L5 | jac | F", 2 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort2_l6_jac_f("MGBench mcsse | q2 | sort 2 | L6 | jac | F", 2 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l4_jac_f("MGBench mcsse | q2 | sort 3 | L4 | jac | F", 2 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l5_jac_f("MGBench mcsse | q2 | sort 3 | L5 | jac | F", 2 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort3_l6_jac_f("MGBench mcsse | q2 | sort 3 | L6 | jac | F", 2 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l4_jac_f("MGBench mcsse | q2 | sort 4 | L4 | jac | F", 2 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l5_jac_f("MGBench mcsse | q2 | sort 4 | L5 | jac | F", 2 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q2_sort4_l6_jac_f("MGBench mcsse | q2 | sort 4 | L6 | jac | F", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l4_jac_f("MGBench mcsse | q1t | sort 0 | L4 | jac | F", 3 , 0, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l5_jac_f("MGBench mcsse | q1t | sort 0 | L5 | jac | F", 3 , 0, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort0_l6_jac_f("MGBench mcsse | q1t | sort 0 | L6 | jac | F", 3 , 0, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l4_jac_f("MGBench mcsse | q1t | sort 1 | L4 | jac | F", 3 , 1, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l5_jac_f("MGBench mcsse | q1t | sort 1 | L5 | jac | F", 3 , 1, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort1_l6_jac_f("MGBench mcsse | q1t | sort 1 | L6 | jac | F", 3 , 1, 6, "jac", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l4_jac_f("MGBench mcsse | q1t | sort 2 | L4 | jac | F", 3 , 2, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l5_jac_f("MGBench mcsse | q1t | sort 2 | L5 | jac | F", 3 , 2, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort2_l6_jac_f("MGBench mcsse | q1t | sort 2 | L6 | jac | F", 3 , 2, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l4_jac_f("MGBench mcsse | q1t | sort 3 | L4 | jac | F", 3 , 3, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l5_jac_f("MGBench mcsse | q1t | sort 3 | L5 | jac | F", 3 , 3, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort3_l6_jac_f("MGBench mcsse | q1t | sort 3 | L6 | jac | F", 3 , 3, 6, "jac", 0.5);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l4_jac_f("MGBench mcsse | q1t | sort 4 | L4 | jac | F", 3 , 4, 4, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l5_jac_f("MGBench mcsse | q1t | sort 4 | L5 | jac | F", 3 , 4, 5, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, DenseVector<double> > mcsse_q1t_sort4_l6_jac_f("MGBench mcsse | q1t | sort 4 | L6 | jac | F", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_spai_f("MGBench mcsse | q1 | sort 0 | L5 | spai | F", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_spai_f("MGBench mcsse | q1 | sort 0 | L6 | spai | F", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_spai_f("MGBench mcsse | q1 | sort 0 | L7 | spai | F", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_spai_f("MGBench mcsse | q1 | sort 1 | L5 | spai | F", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_spai_f("MGBench mcsse | q1 | sort 1 | L6 | spai | F", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_spai_f("MGBench mcsse | q1 | sort 1 | L7 | spai | F", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_spai_f("MGBench mcsse | q1 | sort 2 | L5 | spai | F", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_spai_f("MGBench mcsse | q1 | sort 2 | L6 | spai | F", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_spai_f("MGBench mcsse | q1 | sort 2 | L7 | spai | F", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_spai_f("MGBench mcsse | q1 | sort 3 | L5 | spai | F", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_spai_f("MGBench mcsse | q1 | sort 3 | L6 | spai | F", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_spai_f("MGBench mcsse | q1 | sort 3 | L7 | spai | F", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_spai_f("MGBench mcsse | q1 | sort 4 | L5 | spai | F", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_spai_f("MGBench mcsse | q1 | sort 4 | L6 | spai | F", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_spai_f("MGBench mcsse | q1 | sort 4 | L7 | spai | F", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_spai_f("MGBench mcsse | q2 | sort 0 | L4 | spai | F", 2 , 0, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_spai_f("MGBench mcsse | q2 | sort 0 | L5 | spai | F", 2 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_spai_f("MGBench mcsse | q2 | sort 0 | L6 | spai | F", 2 , 0, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_spai_f("MGBench mcsse | q2 | sort 1 | L4 | spai | F", 2 , 1, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_spai_f("MGBench mcsse | q2 | sort 1 | L5 | spai | F", 2 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_spai_f("MGBench mcsse | q2 | sort 1 | L6 | spai | F", 2 , 1, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_spai_f("MGBench mcsse | q2 | sort 2 | L4 | spai | F", 2 , 2, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_spai_f("MGBench mcsse | q2 | sort 2 | L5 | spai | F", 2 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_spai_f("MGBench mcsse | q2 | sort 2 | L6 | spai | F", 2 , 2, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_spai_f("MGBench mcsse | q2 | sort 3 | L4 | spai | F", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_spai_f("MGBench mcsse | q2 | sort 3 | L5 | spai | F", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_spai_f("MGBench mcsse | q2 | sort 3 | L6 | spai | F", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_spai_f("MGBench mcsse | q2 | sort 4 | L4 | spai | F", 2 , 4, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_spai_f("MGBench mcsse | q2 | sort 4 | L5 | spai | F", 2 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_spai_f("MGBench mcsse | q2 | sort 4 | L6 | spai | F", 2 , 4, 6, "spai_grote", 1.);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_spai_f("MGBench mcsse | q1t | sort 0 | L4 | spai | F", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_spai_f("MGBench mcsse | q1t | sort 0 | L5 | spai | F", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_spai_f("MGBench mcsse | q1t | sort 0 | L6 | spai | F", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_spai_f("MGBench mcsse | q1t | sort 1 | L4 | spai | F", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_spai_f("MGBench mcsse | q1t | sort 1 | L5 | spai | F", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_spai_f("MGBench mcsse | q1t | sort 1 | L6 | spai | F", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_spai_f("MGBench mcsse | q1t | sort 2 | L4 | spai | F", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_spai_f("MGBench mcsse | q1t | sort 2 | L5 | spai | F", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_spai_f("MGBench mcsse | q1t | sort 2 | L6 | spai | F", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_spai_f("MGBench mcsse | q1t | sort 3 | L4 | spai | F", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_spai_f("MGBench mcsse | q1t | sort 3 | L5 | spai | F", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_spai_f("MGBench mcsse | q1t | sort 3 | L6 | spai | F", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_spai_f("MGBench mcsse | q1t | sort 4 | L4 | spai | F", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_spai_f("MGBench mcsse | q1t | sort 4 | L5 | spai | F", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_spai_f("MGBench mcsse | q1t | sort 4 | L6 | spai | F", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l5_sainv_f("MGBench mcsse | q1 | sort 0 | L5 | sainv | F", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l6_sainv_f("MGBench mcsse | q1 | sort 0 | L6 | sainv | F", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort0_l7_sainv_f("MGBench mcsse | q1 | sort 0 | L7 | sainv | F", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l5_sainv_f("MGBench mcsse | q1 | sort 1 | L5 | sainv | F", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l6_sainv_f("MGBench mcsse | q1 | sort 1 | L6 | sainv | F", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort1_l7_sainv_f("MGBench mcsse | q1 | sort 1 | L7 | sainv | F", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l5_sainv_f("MGBench mcsse | q1 | sort 2 | L5 | sainv | F", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l6_sainv_f("MGBench mcsse | q1 | sort 2 | L6 | sainv | F", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort2_l7_sainv_f("MGBench mcsse | q1 | sort 2 | L7 | sainv | F", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l5_sainv_f("MGBench mcsse | q1 | sort 3 | L5 | sainv | F", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l6_sainv_f("MGBench mcsse | q1 | sort 3 | L6 | sainv | F", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort3_l7_sainv_f("MGBench mcsse | q1 | sort 3 | L7 | sainv | F", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l5_sainv_f("MGBench mcsse | q1 | sort 4 | L5 | sainv | F", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l6_sainv_f("MGBench mcsse | q1 | sort 4 | L6 | sainv | F", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1_sort4_l7_sainv_f("MGBench mcsse | q1 | sort 4 | L7 | sainv | F", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l4_sainv_f("MGBench mcsse | q2 | sort 0 | L4 | sainv | F", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l5_sainv_f("MGBench mcsse | q2 | sort 0 | L5 | sainv | F", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort0_l6_sainv_f("MGBench mcsse | q2 | sort 0 | L6 | sainv | F", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l4_sainv_f("MGBench mcsse | q2 | sort 1 | L4 | sainv | F", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l5_sainv_f("MGBench mcsse | q2 | sort 1 | L5 | sainv | F", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort1_l6_sainv_f("MGBench mcsse | q2 | sort 1 | L6 | sainv | F", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l4_sainv_f("MGBench mcsse | q2 | sort 2 | L4 | sainv | F", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l5_sainv_f("MGBench mcsse | q2 | sort 2 | L5 | sainv | F", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort2_l6_sainv_f("MGBench mcsse | q2 | sort 2 | L6 | sainv | F", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l4_sainv_f("MGBench mcsse | q2 | sort 3 | L4 | sainv | F", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l5_sainv_f("MGBench mcsse | q2 | sort 3 | L5 | sainv | F", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort3_l6_sainv_f("MGBench mcsse | q2 | sort 3 | L6 | sainv | F", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l4_sainv_f("MGBench mcsse | q2 | sort 4 | L4 | sainv | F", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l5_sainv_f("MGBench mcsse | q2 | sort 4 | L5 | sainv | F", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q2_sort4_l6_sainv_f("MGBench mcsse | q2 | sort 4 | L6 | sainv | F", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l4_sainv_f("MGBench mcsse | q1t | sort 0 | L4 | sainv | F", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l5_sainv_f("MGBench mcsse | q1t | sort 0 | L5 | sainv | F", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort0_l6_sainv_f("MGBench mcsse | q1t | sort 0 | L6 | sainv | F", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l4_sainv_f("MGBench mcsse | q1t | sort 1 | L4 | sainv | F", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l5_sainv_f("MGBench mcsse | q1t | sort 1 | L5 | sainv | F", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort1_l6_sainv_f("MGBench mcsse | q1t | sort 1 | L6 | sainv | F", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l4_sainv_f("MGBench mcsse | q1t | sort 2 | L4 | sainv | F", 3 , 2, 4, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l5_sainv_f("MGBench mcsse | q1t | sort 2 | L5 | sainv | F", 3 , 2, 5, "sainv", 0.9);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort2_l6_sainv_f("MGBench mcsse | q1t | sort 2 | L6 | sainv | F", 3 , 2, 6, "sainv", 0.9);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l4_sainv_f("MGBench mcsse | q1t | sort 3 | L4 | sainv | F", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l5_sainv_f("MGBench mcsse | q1t | sort 3 | L5 | sainv | F", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort3_l6_sainv_f("MGBench mcsse | q1t | sort 3 | L6 | sainv | F", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l4_sainv_f("MGBench mcsse | q1t | sort 4 | L4 | sainv | F", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l5_sainv_f("MGBench mcsse | q1t | sort 4 | L5 | sainv | F", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > mcsse_q1t_sort4_l6_sainv_f("MGBench mcsse | q1t | sort 4 | L6 | sainv | F", 3 , 4, 6, "sainv", 0.7);
#endif

#if defined HONEI_CUDA && defined HONEI_CUDA_DOUBLE
//JAC
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l5_jac_f("MGBench cuda | q1 | sort 0 | L5 | jac | F", 1 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l6_jac_f("MGBench cuda | q1 | sort 0 | L6 | jac | F", 1 , 0, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort0_l7_jac_f("MGBench cuda | q1 | sort 0 | L7 | jac | F", 1 , 0, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l5_jac_f("MGBench cuda | q1 | sort 1 | L5 | jac | F", 1 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l6_jac_f("MGBench cuda | q1 | sort 1 | L6 | jac | F", 1 , 1, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort1_l7_jac_f("MGBench cuda | q1 | sort 1 | L7 | jac | F", 1 , 1, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort2_l5_jac_f("MGBench cuda | q1 | sort 2 | L5 | jac | F", 1 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort2_l6_jac_f("MGBench cuda | q1 | sort 2 | L6 | jac | F", 1 , 2, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort2_l7_jac_f("MGBench cuda | q1 | sort 2 | L7 | jac | F", 1 , 2, 7, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l5_jac_f("MGBench cuda | q1 | sort 3 | L5 | jac | F", 1 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l6_jac_f("MGBench cuda | q1 | sort 3 | L6 | jac | F", 1 , 3, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort3_l7_jac_f("MGBench cuda | q1 | sort 3 | L7 | jac | F", 1 , 3, 7, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort4_l5_jac_f("MGBench cuda | q1 | sort 4 | L5 | jac | F", 1 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort4_l6_jac_f("MGBench cuda | q1 | sort 4 | L6 | jac | F", 1 , 4, 6, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1_sort4_l7_jac_f("MGBench cuda | q1 | sort 4 | L7 | jac | F", 1 , 4, 7, "jac", 0.5);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort0_l4_jac_f("MGBench cuda | q2 | sort 0 | L4 | jac | F", 2 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort0_l5_jac_f("MGBench cuda | q2 | sort 0 | L5 | jac | F", 2 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort0_l6_jac_f("MGBench cuda | q2 | sort 0 | L6 | jac | F", 2 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort1_l4_jac_f("MGBench cuda | q2 | sort 1 | L4 | jac | F", 2 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort1_l5_jac_f("MGBench cuda | q2 | sort 1 | L5 | jac | F", 2 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort1_l6_jac_f("MGBench cuda | q2 | sort 1 | L6 | jac | F", 2 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort2_l4_jac_f("MGBench cuda | q2 | sort 2 | L4 | jac | F", 2 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort2_l5_jac_f("MGBench cuda | q2 | sort 2 | L5 | jac | F", 2 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort2_l6_jac_f("MGBench cuda | q2 | sort 2 | L6 | jac | F", 2 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort3_l4_jac_f("MGBench cuda | q2 | sort 3 | L4 | jac | F", 2 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort3_l5_jac_f("MGBench cuda | q2 | sort 3 | L5 | jac | F", 2 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort3_l6_jac_f("MGBench cuda | q2 | sort 3 | L6 | jac | F", 2 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort4_l4_jac_f("MGBench cuda | q2 | sort 4 | L4 | jac | F", 2 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort4_l5_jac_f("MGBench cuda | q2 | sort 4 | L5 | jac | F", 2 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q2_sort4_l6_jac_f("MGBench cuda | q2 | sort 4 | L6 | jac | F", 2 , 4, 6, "jac", 0.5);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l4_jac_f("MGBench cuda | q1t | sort 0 | L4 | jac | F", 3 , 0, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l5_jac_f("MGBench cuda | q1t | sort 0 | L5 | jac | F", 3 , 0, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort0_l6_jac_f("MGBench cuda | q1t | sort 0 | L6 | jac | F", 3 , 0, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l4_jac_f("MGBench cuda | q1t | sort 1 | L4 | jac | F", 3 , 1, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l5_jac_f("MGBench cuda | q1t | sort 1 | L5 | jac | F", 3 , 1, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort1_l6_jac_f("MGBench cuda | q1t | sort 1 | L6 | jac | F", 3 , 1, 6, "jac", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l4_jac_f("MGBench cuda | q1t | sort 2 | L4 | jac | F", 3 , 2, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l5_jac_f("MGBench cuda | q1t | sort 2 | L5 | jac | F", 3 , 2, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort2_l6_jac_f("MGBench cuda | q1t | sort 2 | L6 | jac | F", 3 , 2, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l4_jac_f("MGBench cuda | q1t | sort 3 | L4 | jac | F", 3 , 3, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l5_jac_f("MGBench cuda | q1t | sort 3 | L5 | jac | F", 3 , 3, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort3_l6_jac_f("MGBench cuda | q1t | sort 3 | L6 | jac | F", 3 , 3, 6, "jac", 0.5);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l4_jac_f("MGBench cuda | q1t | sort 4 | L4 | jac | F", 3 , 4, 4, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l5_jac_f("MGBench cuda | q1t | sort 4 | L5 | jac | F", 3 , 4, 5, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, DenseVector<double> > cuda_q1t_sort4_l6_jac_f("MGBench cuda | q1t | sort 4 | L6 | jac | F", 3 , 4, 6, "jac", 0.5);

//SPAI
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_spai_f("MGBench cuda | q1 | sort 0 | L5 | spai | F", 1 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_spai_f("MGBench cuda | q1 | sort 0 | L6 | spai | F", 1 , 0, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_spai_f("MGBench cuda | q1 | sort 0 | L7 | spai | F", 1 , 0, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_spai_f("MGBench cuda | q1 | sort 1 | L5 | spai | F", 1 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_spai_f("MGBench cuda | q1 | sort 1 | L6 | spai | F", 1 , 1, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_spai_f("MGBench cuda | q1 | sort 1 | L7 | spai | F", 1 , 1, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_spai_f("MGBench cuda | q1 | sort 2 | L5 | spai | F", 1 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_spai_f("MGBench cuda | q1 | sort 2 | L6 | spai | F", 1 , 2, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_spai_f("MGBench cuda | q1 | sort 2 | L7 | spai | F", 1 , 2, 7, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_spai_f("MGBench cuda | q1 | sort 3 | L5 | spai | F", 1 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_spai_f("MGBench cuda | q1 | sort 3 | L6 | spai | F", 1 , 3, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_spai_f("MGBench cuda | q1 | sort 3 | L7 | spai | F", 1 , 3, 7, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_spai_f("MGBench cuda | q1 | sort 4 | L5 | spai | F", 1 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_spai_f("MGBench cuda | q1 | sort 4 | L6 | spai | F", 1 , 4, 6, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_spai_f("MGBench cuda | q1 | sort 4 | L7 | spai | F", 1 , 4, 7, "spai_grote", 1.);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_spai_f("MGBench cuda | q2 | sort 0 | L4 | spai | F", 2 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_spai_f("MGBench cuda | q2 | sort 0 | L5 | spai | F", 2 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_spai_f("MGBench cuda | q2 | sort 0 | L6 | spai | F", 2 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_spai_f("MGBench cuda | q2 | sort 1 | L4 | spai | F", 2 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_spai_f("MGBench cuda | q2 | sort 1 | L5 | spai | F", 2 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_spai_f("MGBench cuda | q2 | sort 1 | L6 | spai | F", 2 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_spai_f("MGBench cuda | q2 | sort 2 | L4 | spai | F", 2 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_spai_f("MGBench cuda | q2 | sort 2 | L5 | spai | F", 2 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_spai_f("MGBench cuda | q2 | sort 2 | L6 | spai | F", 2 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_spai_f("MGBench cuda | q2 | sort 3 | L4 | spai | F", 2 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_spai_f("MGBench cuda | q2 | sort 3 | L5 | spai | F", 2 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_spai_f("MGBench cuda | q2 | sort 3 | L6 | spai | F", 2 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_spai_f("MGBench cuda | q2 | sort 4 | L4 | spai | F", 2 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_spai_f("MGBench cuda | q2 | sort 4 | L5 | spai | F", 2 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_spai_f("MGBench cuda | q2 | sort 4 | L6 | spai | F", 2 , 4, 6, "spai_grote", 1.);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_spai_f("MGBench cuda | q1t | sort 0 | L4 | spai | F", 3 , 0, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_spai_f("MGBench cuda | q1t | sort 0 | L5 | spai | F", 3 , 0, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_spai_f("MGBench cuda | q1t | sort 0 | L6 | spai | F", 3 , 0, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_spai_f("MGBench cuda | q1t | sort 1 | L4 | spai | F", 3 , 1, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_spai_f("MGBench cuda | q1t | sort 1 | L5 | spai | F", 3 , 1, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_spai_f("MGBench cuda | q1t | sort 1 | L6 | spai | F", 3 , 1, 6, "spai_grote", 1.);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_spai_f("MGBench cuda | q1t | sort 2 | L4 | spai | F", 3 , 2, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_spai_f("MGBench cuda | q1t | sort 2 | L5 | spai | F", 3 , 2, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_spai_f("MGBench cuda | q1t | sort 2 | L6 | spai | F", 3 , 2, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_spai_f("MGBench cuda | q1t | sort 3 | L4 | spai | F", 3 , 3, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_spai_f("MGBench cuda | q1t | sort 3 | L5 | spai | F", 3 , 3, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_spai_f("MGBench cuda | q1t | sort 3 | L6 | spai | F", 3 , 3, 6, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_spai_f("MGBench cuda | q1t | sort 4 | L4 | spai | F", 3 , 4, 4, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_spai_f("MGBench cuda | q1t | sort 4 | L5 | spai | F", 3 , 4, 5, "spai_grote", 1.);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_spai_f("MGBench cuda | q1t | sort 4 | L6 | spai | F", 3 , 4, 6, "spai_grote", 1.);

//SAINV
//q1
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l5_sainv_f("MGBench cuda | q1 | sort 0 | L5 | sainv | F", 1 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l6_sainv_f("MGBench cuda | q1 | sort 0 | L6 | sainv | F", 1 , 0, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort0_l7_sainv_f("MGBench cuda | q1 | sort 0 | L7 | sainv | F", 1 , 0, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l5_sainv_f("MGBench cuda | q1 | sort 1 | L5 | sainv | F", 1 , 1, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l6_sainv_f("MGBench cuda | q1 | sort 1 | L6 | sainv | F", 1 , 1, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort1_l7_sainv_f("MGBench cuda | q1 | sort 1 | L7 | sainv | F", 1 , 1, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l5_sainv_f("MGBench cuda | q1 | sort 2 | L5 | sainv | F", 1 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l6_sainv_f("MGBench cuda | q1 | sort 2 | L6 | sainv | F", 1 , 2, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort2_l7_sainv_f("MGBench cuda | q1 | sort 2 | L7 | sainv | F", 1 , 2, 7, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l5_sainv_f("MGBench cuda | q1 | sort 3 | L5 | sainv | F", 1 , 3, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l6_sainv_f("MGBench cuda | q1 | sort 3 | L6 | sainv | F", 1 , 3, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort3_l7_sainv_f("MGBench cuda | q1 | sort 3 | L7 | sainv | F", 1 , 3, 7, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l5_sainv_f("MGBench cuda | q1 | sort 4 | L5 | sainv | F", 1 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l6_sainv_f("MGBench cuda | q1 | sort 4 | L6 | sainv | F", 1 , 4, 6, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1_sort4_l7_sainv_f("MGBench cuda | q1 | sort 4 | L7 | sainv | F", 1 , 4, 7, "sainv", 0.9);

//q2
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l4_sainv_f("MGBench cuda | q2 | sort 0 | L4 | sainv | F", 2 , 0, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l5_sainv_f("MGBench cuda | q2 | sort 0 | L5 | sainv | F", 2 , 0, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort0_l6_sainv_f("MGBench cuda | q2 | sort 0 | L6 | sainv | F", 2 , 0, 6, "sainv", 0.9);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l4_sainv_f("MGBench cuda | q2 | sort 1 | L4 | sainv | F", 2 , 1, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l5_sainv_f("MGBench cuda | q2 | sort 1 | L5 | sainv | F", 2 , 1, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort1_l6_sainv_f("MGBench cuda | q2 | sort 1 | L6 | sainv | F", 2 , 1, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l4_sainv_f("MGBench cuda | q2 | sort 2 | L4 | sainv | F", 2 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l5_sainv_f("MGBench cuda | q2 | sort 2 | L5 | sainv | F", 2 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort2_l6_sainv_f("MGBench cuda | q2 | sort 2 | L6 | sainv | F", 2 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l4_sainv_f("MGBench cuda | q2 | sort 3 | L4 | sainv | F", 2 , 3, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l5_sainv_f("MGBench cuda | q2 | sort 3 | L5 | sainv | F", 2 , 3, 5, "sainv", 0.9);
//MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort3_l6_sainv_f("MGBench cuda | q2 | sort 3 | L6 | sainv | F", 2 , 3, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l4_sainv_f("MGBench cuda | q2 | sort 4 | L4 | sainv | F", 2 , 4, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l5_sainv_f("MGBench cuda | q2 | sort 4 | L5 | sainv | F", 2 , 4, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q2_sort4_l6_sainv_f("MGBench cuda | q2 | sort 4 | L6 | sainv | F", 2 , 4, 6, "sainv", 0.9);

//q1t
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l4_sainv_f("MGBench cuda | q1t | sort 0 | L4 | sainv | F", 3 , 0, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l5_sainv_f("MGBench cuda | q1t | sort 0 | L5 | sainv | F", 3 , 0, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort0_l6_sainv_f("MGBench cuda | q1t | sort 0 | L6 | sainv | F", 3 , 0, 6, "sainv", 0.7);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l4_sainv_f("MGBench cuda | q1t | sort 1 | L4 | sainv | F", 3 , 1, 4, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l5_sainv_f("MGBench cuda | q1t | sort 1 | L5 | sainv | F", 3 , 1, 5, "sainv", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort1_l6_sainv_f("MGBench cuda | q1t | sort 1 | L6 | sainv | F", 3 , 1, 6, "sainv", 0.5);


MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l4_sainv_f("MGBench cuda | q1t | sort 2 | L4 | sainv | F", 3 , 2, 4, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l5_sainv_f("MGBench cuda | q1t | sort 2 | L5 | sainv | F", 3 , 2, 5, "sainv", 0.9);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort2_l6_sainv_f("MGBench cuda | q1t | sort 2 | L6 | sainv | F", 3 , 2, 6, "sainv", 0.9);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l4_sainv_f("MGBench cuda | q1t | sort 3 | L4 | sainv | F", 3 , 3, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l5_sainv_f("MGBench cuda | q1t | sort 3 | L5 | sainv | F", 3 , 3, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort3_l6_sainv_f("MGBench cuda | q1t | sort 3 | L6 | sainv | F", 3 , 3, 6, "sainv", 0.7);

MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l4_sainv_f("MGBench cuda | q1t | sort 4 | L4 | sainv | F", 3 , 4, 4, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l5_sainv_f("MGBench cuda | q1t | sort 4 | L5 | sainv | F", 3 , 4, 5, "sainv", 0.7);
MGBench<tags::GPU::CUDA, methods::CYCLE::F::V::STATIC, SparseMatrixELL<double> > cuda_q1t_sort4_l6_sainv_f("MGBench cuda | q1t | sort 4 | L6 | sainv | F", 3 , 4, 6, "sainv", 0.7);
#endif
