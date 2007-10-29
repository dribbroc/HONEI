/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBMATH_GUARD_ITERATIVE_REFINEMENT_HH
#define LIBMATH_GUARD_ITERATIVE_REFINEMENT_HH 1

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#include <libmath/conjugate_gradients.hh>
#include <libmath/jacobi.hh>
/**
 * \file
 *
 * Templatized definition and implementation of a Mixed Precision Iterative
 * Refinement solver.
 *
 * \ingroup grpoperations
 */

using namespace std;
namespace honei
{
    /**
     * \brief Solution of LES with Iterative Refinement using CG as inner solver.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct IterativeRefinement
    {
        public:
            template<typename DT1_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                DenseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size(), float(0));
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                defect = Difference<Tag_>::value(defect, right_hand_side);
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(DT1_(1./alpha), defect);
                        }
                        else
                        {
                            Scale<Tag_>::value(DT1_(1./ std::numeric_limits<double>::epsilon()), defect);
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseMatrix<DT1_>::ConstElementIterator i_outer(system_matrix.begin_elements()), i_end(system_matrix.end_elements());
                    typename DenseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                    while(i_outer != i_end)
                    {
                        *i_inner = float(*i_outer);
                        ++i_inner; ++i_outer;
                    }

                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }

                    inner_defect = ConjugateGradients<Tag_>::value(inner_system, inner_defect, eps_inner);
                    //defect = ConjugateGradients<Tag_>::value(system_matrix, defect, eps_inner);

                    typename DenseMatrix<DT1_>::ElementIterator a_outer(system_matrix.begin_elements()), a_end(system_matrix.end_elements());
                    typename DenseMatrix<float>::ConstElementIterator a_inner(inner_system.begin_elements());
                    while(a_outer != a_end)
                    {
                        *a_outer = DT1_(*a_inner);
                        ++a_inner; ++a_outer;
                    }

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = float(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    DenseVector<DT1_> c_scaled = Scale<Tag_>::value(alpha, defect);
                    x_actual = Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    defect = Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }
    };
}

#endif
