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

#include <operator.hh>
#include <operator_list.hh>
#include <honei/util/profiler.hh>
#include <methods.hh>

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
            unsigned long used_iters;
            unsigned long n_pre_smooth;
            unsigned long n_post_smooth;
            unsigned long min_level;
            double eps_relative;
    };

    template<typename MatrixType_, typename VectorType_, typename PreconContType_>
    struct MGUtil
    {
        public:
            MGData<MatrixType_, VectorType_, PreconContType_> load_data(std::string file_base, unsigned long sorting)
            {
                CONTEXT("When creating MGData from file:");

                std::vector<MatrixType_> system;
                std::vector<PreconContType_> P;
                std::vector<VectorType_> b;
                std::vector<VectorType_> x;
                std::vector<VectorType_> d;
                std::vector<VectorType_> c;
                std::vector<VectorType_> temp_0;
                std::vector<VectorType_> temp_1;

                //TODO: load and push back data

                MGData<MatrixType_, VectorType_, PreconContType_> result(system, P, b, x, d, c, temp_0, temp_1);
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
                    cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_>(
                                data.A.at(data.min_level),
                                b.at(data.min_level),
                                x.at(data.min_level),
                                data.max_iters_coarse,
                                data.used_iters,
                                data.eps_relative) );
                }
                else
                {
                    ///Presmoothing
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level),
                                data.P.at(level),
                                b.at(level),
                                x.at(level),
                                data.temp_0.at(level),
                                data.temp_1.at(level),
                                data.n_pre_smooth) );

                    ///Defect
                    cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.d.at(level), b.at(level), data.A.at(level), x.at(level)));

                    ///Restriction
                    cycle.push_back(new TransferOperator<ResType_, MatrixType_, VectorType_>(data.d.at(level - 1) , data.d.at(level), data.resmat.at(level)));

                    ///Recursion
                    ///all vectors in c have to be initialised with 0
                    _coarse_correction(data.c, data.d, level - 1, cycle, data);

                    ///Prolongation
                    cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(data.c.at(level + 1) , data.c.at(level), data.prolmat.at(level)));
                    cycle.push_back(new SumOperator<Tag_, VectorType_>(data.x.at(level + 1), data.c.at(level + 1)));

                    ///Postsmoothing
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level),
                                data.P.at(level),
                                data.b.at(level),
                                data.x.at(level),
                                data.temp_0.at(level),
                                data.temp_1.at(level),
                                data.n_post_smooth) );
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
                    cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_>(
                                data.A.at(data.min_level),
                                b.at(data.min_level),
                                x.at(data.min_level),
                                data.max_iters_coarse,
                                data.used_iters,
                                data.eps_relative) );
                }
                else
                {
                    ///Presmoothing
                    cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                data.A.at(level),
                                data.P.at(level),
                                b.at(level),
                                x.at(level),
                                data.temp_0.at(level),
                                data.temp_1.at(level),
                                data.n_pre_smooth) );

                    ///Restriction
                    cycle.push_back(new TransferOperator<ResType_, MatrixType_, VectorType_>(data.d.at(level - 1) , data.d.at(level), data.resmat.at(level)));

                    ///Recursion
                    ///all vectors in c have to be initialised with 0
                    _coarse_correction(x, b, level - 1, cycle, data);

                    ///Prolongation
                    cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(x.at(level + 1) , x.at(level), data.prolmat.at(level)));
                    cycle.push_back(new SumOperator<Tag_, VectorType_>(data.x.at(level + 1), x.at(level + 1)));

                    ///Postsmoothing
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
                _build_cycle(data.x, data.b, data.A.size() - 1, cycle, data);
                return cycle;
            }
    };
}

#endif
