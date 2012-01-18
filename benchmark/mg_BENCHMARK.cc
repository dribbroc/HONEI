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
            file += "/honei/math/testdata/poisson_advanced2/";

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


            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, PreconContType_, double >  data(MGUtil<Tag_,
                    SparseMatrixELL<double>,
                    DenseVector<double>,
                    SparseMatrixELL<double>,
                    PreconContType_,
                    MatrixIO<io_formats::ELL>,
                    VectorIO<io_formats::EXP>,
                    double>::load_data(file, _levels, _damping, _precon));
            MGUtil<Tag_,
                SparseMatrixELL<double>,
                DenseVector<double>,
                SparseMatrixELL<double>,
                PreconContType_,
                MatrixIO<io_formats::ELL>,
                VectorIO<io_formats::EXP>,
                double>::configure(data, 100, 10, 4, 4, 1, double(1e-8));

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

//JAC
//q1
MGBench<tags::CPU::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > sse_q1_sort2_l7_jac_v("MGBench mcsse | q1 | sort 2 | L7 | jac | V", 1 , 2, 7, "jac", 0.7);
MGBench<tags::CPU::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > sse_q1_sort2_l8_jac_v("MGBench mcsse | q1 | sort 2 | L8 | jac | V", 1 , 2, 8, "jac", 0.7);
MGBench<tags::CPU::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > sse_q1_sort2_l9_jac_v("MGBench mcsse | q1 | sort 2 | L9 | jac | V", 1 , 2, 9, "jac", 0.7);


MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l7_jac_v("MGBench mcsse | q1 | sort 2 | L7 | jac | V", 1 , 2, 7, "jac", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l8_jac_v("MGBench mcsse | q1 | sort 2 | L8 | jac | V", 1 , 2, 8, "jac", 0.7);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, DenseVector<double> > mcsse_q1_sort2_l9_jac_v("MGBench mcsse | q1 | sort 2 | L9 | jac | V", 1 , 2, 9, "jac", 0.7);

