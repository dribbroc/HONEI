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
#define BENCHMARK_INSTRUMENTATION 1
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


template <typename Tag_, typename CycleType_, typename MT_, typename PreconContType_>
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
            //file += "/honei/math/testdata/poisson_advanced4/";
            file += "/honei/math/testdata/bench1-fies/";

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
            bool adaptive = false;

            _damping = adaptive ? double(1) : _damping;

            MGData<MT_, DenseVector<double>, MT_, PreconContType_, double >  data(MGUtil<Tag_,
                    MT_,
                    DenseVector<double>,
                    MT_,
                    PreconContType_,
                    MatrixIO<io_formats::ELL>,
                    VectorIO<io_formats::EXP>,
                    double>::load_data(file, _levels, _damping, _precon));

            MGUtil<Tag_,
                MT_,
                DenseVector<double>,
                MT_,
                PreconContType_,
                MatrixIO<io_formats::ELL>,
                VectorIO<io_formats::EXP>,
                double>::configure(data, 2000, 10, 4, 4, 1, double(1e-8));

            OperatorList ol(
                    MGCycleCreation<Tag_,
                    CycleType_,
                    //SuperLU,
                    BiCGStabSolver<Tag_, methods::VAR>,
                    //CG<Tag_, methods::NONE>,
                    RISmoother<Tag_>,
                    //BiCGStabSmoother<Tag_>,
                    Restriction<Tag_, methods::PROLMAT>,
                    Prolongation<Tag_, methods::PROLMAT>,
                    double>::value(data)
                    );


            BENCHRESET();

            BENCHMARK(
                    (MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data, ol));
#ifdef HONEI_CUDA
                    if (Tag_::tag_value == tags::tv_gpu_cuda)
                    cuda::GPUPool::instance()->flush();
#endif
                    );

            evaluate(globalBenchmarkInfo);

            std::cout << data.used_iters << " iterations used." << std::endl;

            std::cout<<"FLOP: "<<globalBenchmarkInfo.flops<<" Load: "<<globalBenchmarkInfo.load<<" Store: "<<globalBenchmarkInfo.store<<std::endl;
        }
};

MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixCSR<double>, DenseVector<double> > mcsse_q1_sort2_l8_spai_jac_v("MGBench mcsse | q1 | sort 0 | L8 | jac | V", 1 , 0, 8, "jac", 0.5);
MGBench<tags::CPU::MultiCore::SSE, methods::CYCLE::V::STATIC, SparseMatrixCSR<double>, SparseMatrixCSR<double> > mcsse_q1_sort2_l8_spai_grote_v("MGBench mcsse | q1 | sort 0 | L8 | spai_grote | V", 1 , 0, 8, "spai_grote", 1.);

MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double>, DenseVector<double> > cuda_q1_sort2_l8_spai_jac_v("MGBench cuda | q1 | sort 0 | L8 | jac | V", 1 , 0, 8, "jac", 0.5);
MGBench<tags::GPU::CUDA, methods::CYCLE::V::STATIC, SparseMatrixELL<double>, SparseMatrixELL<double> > cuda_q1_sort2_l8_spai_grote_v("MGBench cuda | q1 | sort 0 | L8 | spai_grote | V", 1 , 0, 8, "spai_grote", 1.);

//MGBench<tags::CPU::Generic, methods::CYCLE::V::STATIC, SparseMatrixELL<double>, DenseVector<double> > generic_q1_sort2_l8_spai_jac_v("MGBench generic | q1 | sort 0 | L8 | jac | V", 1 , 0, 8, "jac", 0.5);
//MGBench<tags::CPU::Generic, methods::CYCLE::V::STATIC, SparseMatrixELL<double>, SparseMatrixELL<double> > generic_q1_sort2_l8_spai_grote_v("MGBench generic | q1 | sort 0 | L8 | spai_grote | V", 1 , 0, 8, "spai_grote", 1.);

