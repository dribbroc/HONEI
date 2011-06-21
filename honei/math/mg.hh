/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
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


#ifndef MATH_GUARD_MG_HH
#define MATH_GUARD_MG_HH 1

#include <honei/math/operator.hh>
#include <honei/math/operator_list.hh>
#include <honei/util/profiler.hh>
#include <honei/math/methods.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/transposition.hh>
#include <honei/math/vector_io.hh>

namespace honei
{
    template<typename Tag_, typename NormType_>
    struct MG
    {
        public:
            template<typename MatrixType_, typename VectorType_, typename DT_>
            static void value(MatrixType_ & A, VectorType_ & b, VectorType_ & x, OperatorList & cycle, unsigned long max_iters, unsigned long & used_iters, DT_ eps_relative)
            {
                CONTEXT("When solving linear system with MG :");
                ASSERT(cycle.size() > 0, "OperatorList is empty!");
                PROFILER_START("MGSolver");

                VectorType_ r(x.size());
                Defect<Tag_>::value(r, b, A, x);
                DT_ rnorm_initial(NormType_::value(r));
                DT_ rnorm_current(1e16);

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                    Defect<Tag_>::value(r, b, A, x);
                    rnorm_current = NormType_::value(r);

                    used_iters = i + 1;

                    if(rnorm_current < eps_relative * rnorm_initial)
                        break;
                }
                PROFILER_STOP("MGSolver");
            }
    };

    struct MGSmoother
    {
        public:
            static void value(OperatorList & cycle, unsigned long max_iters)
            {
                CONTEXT("When smoothing linear system with MG :");
                ASSERT(cycle.size() > 0, "OperatorList is empty!");
                PROFILER_START("MGSmoother");

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                }
                PROFILER_STOP("MGSmoother");
            }
    };

    template<typename MatrixType_, typename VectorType_, typename PreconContType_>
    struct MGData
    {
        public:
            MGData(std::vector<MatrixType_> & systems,
                   std::vector<MatrixType_> & resmats,
                   std::vector<MatrixType_> & prolmats,
                   std::vector<PreconContType_> & precons,
                   std::vector<VectorType_> & rhss,
                   std::vector<VectorType_> & xs,
                   std::vector<VectorType_> & ds,
                   std::vector<VectorType_> & cs,
                   std::vector<VectorType_> & temp0s,
                   std::vector<VectorType_> & temp1s,
                   unsigned long p_max_iters_coarse,
                   unsigned long & p_used_iters,
                   unsigned long p_n_pre_smooth,
                   unsigned long p_n_post_smooth,
                   unsigned long p_min_level,
                   double p_eps_rel
                   ) :
                A(systems),
                resmat(resmats),
                prolmat(prolmats),
                P(precons),
                b(rhss),
                x(xs),
                d(ds),
                c(cs),
                temp_0(temp0s),
                temp_1(temp1s),
                max_iters_coarse(p_max_iters_coarse),
                used_iters(p_used_iters),
                n_pre_smooth(p_n_pre_smooth),
                n_post_smooth(p_n_post_smooth),
                min_level(p_min_level),
                eps_relative(p_eps_rel)
            {
                CONTEXT("When creating MGData:");
                ASSERT(p_min_level <= systems.size(), "Minimum level is larger then number of levels!");
            }

            std::vector<MatrixType_> A;
            std::vector<MatrixType_> resmat;
            std::vector<MatrixType_> prolmat;
            std::vector<PreconContType_> P;
            std::vector<VectorType_> b;
            std::vector<VectorType_> x;
            std::vector<VectorType_> d;
            std::vector<VectorType_> c;
            std::vector<VectorType_> temp_0;
            std::vector<VectorType_> temp_1;
            unsigned long max_iters_coarse;
            unsigned long & used_iters;
            unsigned long n_pre_smooth;
            unsigned long n_post_smooth;
            unsigned long min_level;
            double eps_relative;
    };

    template<typename Tag_, typename MatrixType_, typename VectorType_, typename PreconContType_, typename MatIOType_, typename VecIOType_, typename DT_>
    struct MGUtil
    {
        public:
            static MGData<MatrixType_, VectorType_, PreconContType_> load_data(std::string file_base, unsigned long max_level)
            {
                CONTEXT("When creating MGData from file(s):");

                std::vector<MatrixType_> A;
                std::vector<MatrixType_> Prol;
                std::vector<MatrixType_> Res;
                std::vector<PreconContType_> P;
                std::vector<VectorType_> b;
                std::vector<VectorType_> x;
                std::vector<VectorType_> d;
                std::vector<VectorType_> c;
                std::vector<VectorType_> temp_0;
                std::vector<VectorType_> temp_1;

                //TODO: load and push back data

                std::string A_name(file_base);
                std::string Prol_name(file_base);
                std::string b_name(file_base);
                std::string x_name(file_base);

                A_name += "A_";
                Prol_name += "prol_";
                b_name += "rhs_";
                x_name += "init_";

                for(unsigned long i(0) ; i < max_level ; ++i)
                {
                    ///get system-matrix for level i+1
                    std::string local_A_name(A_name);
                    local_A_name += stringify(i + 1);
                    local_A_name += ".ell";
                    MatrixType_ local_A(MatrixIO<MatIOType_>::read_matrix(local_A_name, DT_(0)));
                    A.push_back(local_A);

                    ///get Prolmat P_{i+1}^{i+2}
                    if(i == 0)
                    {
                        SparseMatrix<DT_> local_preProl(1,1);
                        MatrixType_ local_Prol(local_preProl);
                        Prol.push_back(local_Prol);
                    }
                    else
                    {
                        std::string local_Prol_name(Prol_name);
                        local_Prol_name += stringify(i + 1);
                        local_Prol_name += ".ell";
                        MatrixType_ local_Prol(MatrixIO<MatIOType_>::read_matrix(local_Prol_name, DT_(0)));
                        Prol.push_back(local_Prol);
                    }

                    ///get Resmat R_{i+2}^{i+1} = (P_{i+1}^{i+2})^T
                    SparseMatrix<DT_> local_preProl(Prol.at(i));
                    SparseMatrix<DT_> local_preRes(Prol.at(i).columns(), Prol.at(i).rows());
                    Transposition<Tag_>::value(local_preProl, local_preRes);
                    MatrixType_ local_Res(local_preRes);
                    Res.push_back(local_Res);

                    ///get vectors for level max_level
                    if(i == max_level - 1)
                    {
                        std::string local_b_name(b_name);
                        local_b_name += stringify(i + 1);
                        VectorType_ max_b(VectorIO<VecIOType_>::read_vector(local_b_name, DT_(0)));
                        b.push_back(max_b);

                        std::string local_x_name(x_name);
                        local_x_name += stringify(i + 1);
                        VectorType_ max_x(VectorIO<VecIOType_>::read_vector(local_x_name, DT_(0)));
                        x.push_back(max_x);
                    }
                    else
                    {
                        ///get vectors for level i+1
                        VectorType_ zero(A.at(i).rows(), DT_(0));
                        b.push_back(zero.copy());
                        x.push_back(zero.copy());
                        d.push_back(zero.copy());
                        c.push_back(zero.copy());
                        temp_0.push_back(zero.copy());
                        temp_1.push_back(zero.copy());
                    }

                    //dirtymost preconhack
                    if(typeid(PreconContType_) != typeid(methods::NONE))
                    {
                        std::string P_name(file_base);
                        P_name += "A_spai_"; //TODO: we should rename all precalculated precons to "P_"; what if there are more types?
                        //TODO
                    }
                }

                unsigned long used_iters(0);
                MGData<MatrixType_, VectorType_, PreconContType_> result(A, Prol, Res, P, b, x, d, c, temp_0, temp_1, 0, used_iters, 0, 0, 0, double(0.));
                return result;
            }
    };

    //TODO: think of more flexible way (more than one smoothertype, ...) ; later, assembly is done here!
    template<typename Tag_,
             typename CycleShape_,
             typename CoarseGridSolverType_,
             typename SmootherType_,
             typename ResType_,
             typename ProlType_,
             typename DT_>
    struct MGCycleProcessing
    {
    };

    //specialise by cycle-shape
    template<typename Tag_,
             typename CoarseGridSolverType_,
             typename SmootherType_,
             typename ResType_,
             typename ProlType_,
             typename DT_>
    struct MGCycleProcessing<Tag_,
                             methods::CYCLE::V::STATIC,
                             CoarseGridSolverType_,
                             SmootherType_,
                             ResType_,
                             ProlType_,
                             DT_>
    {
        private:
            template<typename MatrixType_, typename VectorType_, typename PreconContType_>
            static void _build_cycle(
                    std::vector<VectorType_> & x,
                    std::vector<VectorType_> & b,
                    unsigned long level,
                    OperatorList & cycle,
                    MGData<MatrixType_, VectorType_, PreconContType_> & data)
            {
                if(level == data.min_level)
                {
                    std::cout << "Solver Accessing " << data.min_level - 1 << std::endl;
                    cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_>(
                                data.A.at(data.min_level - 1),
                                b.at(data.min_level - 1),
                                x.at(data.min_level - 1),
                                data.max_iters_coarse,
                                data.used_iters,
                                data.eps_relative) );
                }
                else
                {
                    ///Presmoothing
                    std::cout << "Smoother Accessing " << level - 1 << std::endl;
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level - 1),
                                data.P.at(level - 1),
                                b.at(level - 1),
                                x.at(level - 1),
                                data.temp_0.at(level - 1),
                                data.temp_1.at(level - 1),
                                data.n_pre_smooth) );

                    ///Defect
                    std::cout << "Defect Accessing " << level - 1 << std::endl;
                    cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.d.at(level - 1), b.at(level - 1), data.A.at(level - 1), x.at(level - 1)));

                    ///Restriction
                    std::cout << " Restrict Accessing " << level - 2 << std::endl;
                    std::cout << " Restrict Accessing " << level - 1 << std::endl;
                    cycle.push_back(new TransferOperator<ResType_, MatrixType_, VectorType_>(data.d.at(level - 2) , data.d.at(level - 1), data.resmat.at(level - 1)));

                    ///Recursion
                    ///all vectors in c have to be initialised with 0
                    _coarse_correction(data.c, data.d, level - 1, cycle, data);

                    ///Prolongation
                    /*std::cout << "Prol Accessing " << level << std::endl;
                    std::cout << "Prol Accessing " << level - 1<< std::endl;
                    cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(data.c.at(level) , data.c.at(level - 1), data.prolmat.at(level - 1)));

                    std::cout << "Sum Accessing " << level << std::endl;
                    cycle.push_back(new SumOperator<Tag_, VectorType_>(data.x.at(level), data.c.at(level)));

                    ///Postsmoothing
                    std::cout << "Smoother Accessing " << level - 1 << std::endl;
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level - 1),
                                data.P.at(level - 1),
                                data.b.at(level - 1),
                                data.x.at(level - 1),
                                data.temp_0.at(level - 1),
                                data.temp_1.at(level - 1),
                                data.n_post_smooth) );*/
                }
            }

            template<typename MatrixType_, typename VectorType_, typename PreconContType_>
            static void _coarse_correction(
                    std::vector<VectorType_> & x,
                    std::vector<VectorType_> & b,
                    unsigned long level,
                    OperatorList & cycle,
                    MGData<MatrixType_, VectorType_, PreconContType_> & data)
            {
                if(level == data.min_level)
                {
                    std::cout << "Solver Accessing " << level - 1 << std::endl;
                    cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_>(
                                data.A.at(data.min_level - 1),
                                b.at(data.min_level - 1),
                                x.at(data.min_level - 1),
                                data.max_iters_coarse,
                                data.used_iters,
                                data.eps_relative) );

                    ///Prolongation
                    std::cout << "Prol Accessing " << data.min_level << std::endl;
                    std::cout << "Prol Accessing " << data.min_level - 1 << std::endl;
                    cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(x.at(data.min_level) , x.at(data.min_level - 1), data.prolmat.at(data.min_level - 1)));
                    std::cout << "Sum Accessing " << data.min_level << std::endl;
                    cycle.push_back(new SumOperator<Tag_, VectorType_>(data.x.at(data.min_level), x.at(data.min_level)));

                    ///Postsmoothing
                    std::cout << "Smoother Accessing " << data.min_level << std::endl;
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(data.min_level),
                                data.P.at(data.min_level),
                                b.at(data.min_level),
                                x.at(data.min_level),
                                data.temp_0.at(data.min_level),
                                data.temp_1.at(data.min_level),
                                data.n_post_smooth) );
                }
                else
                {
                    ///Presmoothing
                    std::cout << "Smoother Accessing " << level - 1 << std::endl;
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level - 1),
                                data.P.at(level - 1),
                                b.at(level - 1),
                                x.at(level - 1),
                                data.temp_0.at(level - 1),
                                data.temp_1.at(level - 1),
                                data.n_pre_smooth) );

                    ///Restriction
                    std::cout << "Res Accessing " << level - 2 << std::endl;
                    std::cout << "Res Accessing " << level - 1 << std::endl;
                    cycle.push_back(new TransferOperator<ResType_, MatrixType_, VectorType_>(data.d.at(level - 2) , data.d.at(level - 1), data.resmat.at(level - 1)));

                    ///Recursion
                    ///all vectors in c have to be initialised with 0
                    _coarse_correction(x, b, level - 1, cycle, data);

                    ///Prolongation
                    std::cout << "Prol Accessing " << level << std::endl;
                    std::cout << "Prol Accessing " << level - 1 << std::endl;
                    cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(x.at(level) , x.at(level - 1), data.prolmat.at(level - 1)));
                    std::cout << "Sum Accessing " << level << std::endl;
                    cycle.push_back(new SumOperator<Tag_, VectorType_>(data.x.at(level), x.at(level)));

                    ///Postsmoothing
                    std::cout << "Smoother Accessing " << level << std::endl;
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level),
                                data.P.at(level),
                                b.at(level),
                                x.at(level),
                                data.temp_0.at(level),
                                data.temp_1.at(level),
                                data.n_post_smooth) );
                }
            }

        public:
            template<typename MatrixType_, typename VectorType_, typename PreconContType_>
            static OperatorList value(MGData<MatrixType_, VectorType_, PreconContType_> & data)
            {
                CONTEXT("When evaluating MGCycleProcessing:");

                OperatorList cycle;
                _build_cycle(data.x, data.b, data.A.size(), cycle, data);
                return cycle;
            }
    };
}

#endif
