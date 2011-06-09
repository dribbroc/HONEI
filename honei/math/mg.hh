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

namespace honei
{
    template<typename Tag_, typename NormType_>
    struct MG
    {
        public:
            template<typename MatrixType_, typename VectorType_, typename DT_>
            static void value(MatrixType_ & A, VectorType_ & b, VectorType_ & x, OperatorList & cycle, unsigned long max_iters, unsigned long & used_iters, DT_ eps_relative)
            {
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
            }
    };

    struct MGSmoother
    {
        public:
            static void value(OperatorList & cycle, unsigned long max_iters)
            {
                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                }
            }
    };


    //TODO: think of more flexible way (more than one smoothertype, ...) ; later, assembly is done here!
    template<typename Tag_,
             typename CycleShape_,
             typename NormType_,
             typename CoarseGridSolverType_,
             typename SmootherType_,
             typename PreconContType_,
             typename ResType_,
             typename ProlType_,
             typename DT_>
    struct MGCycleDescriptor
    {
    };

    //specialise by cycle-shape
    template<typename Tag_,
             typename NormType_,
             typename CoarseGridSolverType_,
             typename SmootherType_,
             typename PreconContType_,
             typename ResType_,
             typename ProlType_
             typename DT_>
    struct MGCycleDescriptor<Tag_,
                             cycle::V::STATIC,
                             NormType_,
                             CoarseGridSolverType_,
                             SmootherType_,
                             ResType_,
                             ProlType_,
                             DT_>
    {
        static void value(std::string filename,
                          std::vector<MatrixType_> & A_level,
                          std::vector<MatrixType_> & prolmat_level,
                          std::vector<MatrixType_> & resmat_level,
                          std::vector<PreconContType_> & precon_level,
                          std::vector<VectorType_> & b_level,
                          std::vector<VectorType_> & x_level,
                          std::vector<VectorType_> & c_level,
                          std::vector<VectorType_> & d_level,
                          std::vector<VectorType_> & temp_0_level,
                          std::vector<VectorType_> & temp_1_level,
                          OperatorList & cycle,
                          unsigned long n_pre_smooth,
                          unsigned long n_post_smooth,
                          unsigned long & used_iters,
                          unsigned long min_level = 1,
                          DT_ eps_relative = 1e-8
                          )
        {
            //load data
            MGUtil<MatrixType_, VectorType_, PreconContType_>::load_data(filename,
                                                                         A_level,
                                                                         prolmat_level,
                                                                         resmat_level,
                                                                         precon_level,
                                                                         b_level,
                                                                         x_level,
                                                                         c_level,
                                                                         temp_0_level,
                                                                         temp_1_level);

            //build up operator list
            unsigned long max_level = (unsigned long)A_level.size();

            ///DESCENT
            ///first PRESMOOTHING operator
            cycle.push_back(new SmootherOperator<SmootherType_,
                                                 MatrixType_,
                                                 VectorType_,
                                                 PreconContType_>(A_level[max_level],
                                                                  precon_level[max_level],
                                                                  d_level[max_level],
                                                                  c_level[max_level],
                                                                  temp_0_level[max_level],
                                                                  temp_1_level[max_level],
                                                                  n_pre_smooth) );
            ///SUM operator (add correction to solution)
            cycle.push_back(new SumOperator<Tag_, VectorType_>(x_level[max_level], c_level[max_level]));

            ///DEFECT operator
            cycle.push_back(DefectOperator<Tag_, MatrixType_, VectorType_>(d_level[i], b_level[i], A_level[i], x_level[i]));

            for(unsigned long i(max_level - 1) ; i > min_level ; --i)
            {
                ///PRESMOOTHING operator
                cycle.push_back(new SmootherOperator<SmootherType_,
                                                     MatrixType_,
                                                     VectorType_,
                                                     PreconContType_>(A_level[i],
                                                                      precon_level[i],
                                                                      d_level[i],
                                                                      x_level[i],
                                                                      temp_0_level[i],
                                                                      temp_1_level[i],
                                                                      n_pre_smooth) );
                ///DEFECT operator
                cycle.push_back(DefectOperator<Tag_, MatrixType_, VectorType_>(d_level[i], b_level[i], A_level[i], x_level[i]));

                ///RESTRICTION operator
                cycle.push_back(new TransferOperator<ResType_, MatrixType_, VectorType_>(d_level[i - 1] , d_level[i], resmat[i]));
                //-> d is equivalent to rhs on all levels => we dont need an b-stdvector
            }

            ///COARSE CORRECTION
            ///COARSE GRID SOLVER operator
            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_>(A_level[min_level],
                                                                                                d_level[min_level],
                                                                                                x_level[min_level],
                                                                                                max_iters,
                                                                                                used_iters,
                                                                                                eps_relative) );

            ///RISE
            for(unsigned long i(min_level) ; i < max_level ; ++i)
            {
                ///PROLONGATION operator
                cycle.push_back(new TransferOperator<ProlType_, MatrixType_, VectorType_>(c_level[i + 1] , x_level[i], prolmat[i]));

                ///SUM operator (add correction TODO: adaptive correction)
                cycle.push_back(new SumOperator<Tag_, VectorType_>(x_level[max_level], c_level[max_level]));


                ///POSTMOOTHING operator
                cycle.push_back(new SmootherOperator<SmootherType_,
                                                     MatrixType_,
                                                     VectorType_,
                                                     PreconContType_>(A_level[i],
                                                                      precon_level[i],
                                                                      d_level[i],
                                                                      x_level[i],
                                                                      temp_0_level[i],
                                                                      temp_1_level[i],
                                                                      n_pre_smooth) );
                ///DEFECT operator
                cycle.push_back(DefectOperator<Tag_, MatrixType_, VectorType_>(temp_0_level[i], d_level[i], A_level[i], x_level[i]));

                //BUT NOW?

            }


        }
    };

}

#endif
