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

#include <honei/math/conjugate_gradients.hh>
#include <honei/math/jacobi.hh>
#include <honei/math/methods.hh>

using namespace methods;
namespace honei
{
    template<typename Method_, typename Tag_>
    struct IterativeRefinement
    {
    };

    /**
     * \brief Solution of LES with Iterative Refinement using CG as inner solver.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object (the solution) is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct IterativeRefinement<methods::CG, Tag_>
    {
        public:
            template<typename DT1_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                CONTEXT("When solving dense LES with iterative refinement (CG):");

                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                DenseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size(), float(0));
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                Difference<Tag_>::value(defect, right_hand_side);

                Scale<Tag_>::value(defect, DT1_(-1.));
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Do conversion of system matrix once
                typename DenseMatrix<DT1_>::ConstElementIterator i_outer(system_matrix.begin_elements()), i_end(system_matrix.end_elements());
                typename DenseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                while(i_outer != i_end)
                {
                    *i_inner = float(*i_outer);
                    ++i_inner; ++i_outer;
                }


                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(defect, DT1_(1./alpha));
                        }
                        else
                        {
                            Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }

                    inner_defect = ConjugateGradients<Tag_, methods::NONE>::value(inner_system, inner_defect, eps_inner);

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = DT1_(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    Scale<Tag_>::value(defect, alpha);
                    Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }

            ///Banded System:
            template<typename DT1_>
                static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
                {
                    CONTEXT("When solving banded LES with iterative refinement (CG):");

                    DT1_ alpha = DT1_(1.0);
                    DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                    //inner allocation:
                    BandedMatrix<float> inner_system(right_hand_side.size());
                    DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                    ///Initialize defect and its norm
                    //TODO: use libla/residual(?)

                    DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    Scale<Tag_>::value(defect, DT1_(-1.));
                    DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                    unsigned long iter_number = 0;
                    ///Do conversion of system matrix once:
                    convert(inner_system, system_matrix);

                    ///Main loop:
                    do
                    {
                        ///Scale defect:
                        if(iter_number != 0)
                        {
                            if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                            {
                                Scale<Tag_>::value(defect, DT1_(1./alpha));
                            }
                            else
                            {
                                Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                            }
                        }

                        ///Do conversion and solve inner system:
                        typename DenseVector<DT1_>::ConstElementIterator d_outer(defect.begin_elements()), d_end(defect.end_elements());
                        typename DenseVector<float>::ElementIterator d_inner(inner_defect.begin_elements());
                        while(d_outer != d_end )
                        {
                            *d_inner = float(*d_outer);
                            ++d_inner; ++d_outer;
                        }

                        inner_defect = ConjugateGradients<Tag_, methods::NONE>::value(inner_system, inner_defect, eps_inner);

                        //reconvert
                        convert(defect, inner_defect);

                        ///Update solution:
                        Scale<Tag_>::value(defect, alpha);
                        Sum<Tag_>::value(x_actual, defect);

                        defect = Product<Tag_>::value(system_matrix, x_actual);
                        Difference<Tag_>::value(defect, right_hand_side);

                        alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                        ++iter_number;



                    }
                    while(alpha < eps_outer*initial_defectnorm);

                    return x_actual;
                }

            ///sparse system
            template<typename DT1_>
            static DenseVector<DT1_> value(SparseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                CONTEXT("When solving sparse LES with iterative refinement (CG):");

                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                SparseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size());
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                Difference<Tag_>::value(defect, right_hand_side);

                Scale<Tag_>::value(defect, DT1_(-1.));
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Do conversion of system matrix once
                typename SparseMatrix<DT1_>::NonZeroElementIterator i_outer(system_matrix.begin_non_zero_elements()), i_end(system_matrix.end_non_zero_elements());
                typename SparseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                while(i_outer != i_end)
                {
                    while(i_inner.index() != i_outer.index())
                    {
                        ++i_inner;
                    }
                    *i_inner = float(*i_outer);
                    ++i_outer;
                }

                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(defect, DT1_(1./alpha));
                        }
                        else
                        {
                            Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }
                    DenseVector<float> defect_result(inner_defect.copy());
                    ConjugateGradients<Tag_, methods::NONE>::value(inner_system, inner_defect, defect_result, eps_inner);
                    inner_defect = defect_result;

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = DT1_(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    Scale<Tag_>::value(defect, alpha);
                    Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }


    };
    /**
     * \brief Solution of LES with Iterative Refinement using PCG- Jacobi as inner solver.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object (the solution) is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct IterativeRefinement<methods::PCG::JAC, Tag_>
    {
        public:
            template<typename DT1_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                CONTEXT("When solving dense LES with iterative refinement (CG):");

                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                DenseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size(), float(0));
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                Difference<Tag_>::value(defect, right_hand_side);

                Scale<Tag_>::value(defect, DT1_(-1.));
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Do conversion of system matrix once
                typename DenseMatrix<DT1_>::ConstElementIterator i_outer(system_matrix.begin_elements()), i_end(system_matrix.end_elements());
                typename DenseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                while(i_outer != i_end)
                {
                    *i_inner = float(*i_outer);
                    ++i_inner; ++i_outer;
                }


                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(defect, DT1_(1./alpha));
                        }
                        else
                        {
                            Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }
                    DenseVector<float> defect_result(inner_defect.copy());
                    ConjugateGradients<Tag_, methods::JAC>::value(inner_system, inner_defect, defect_result, eps_inner);
                    inner_defect = defect_result;

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = DT1_(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    Scale<Tag_>::value(defect, alpha);
                    Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }

            ///Banded System:
            template<typename DT1_>
                static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
                {
                    CONTEXT("When solving banded LES with iterative refinement (CG):");

                    DT1_ alpha = DT1_(1.0);
                    DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                    //inner allocation:
                    BandedMatrix<float> inner_system(right_hand_side.size());
                    DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                    ///Initialize defect and its norm
                    //TODO: use libla/residual(?)

                    DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    Scale<Tag_>::value(defect, DT1_(-1.));
                    DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                    unsigned long iter_number = 0;
                    ///Do conversion of system matrix once:
                    convert(inner_system, system_matrix);

                    ///Main loop:
                    do
                    {
                        ///Scale defect:
                        if(iter_number != 0)
                        {
                            if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                            {
                                Scale<Tag_>::value(defect, DT1_(1./alpha));
                            }
                            else
                            {
                                Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                            }
                        }

                        ///Do conversion and solve inner system:
                        typename DenseVector<DT1_>::ConstElementIterator d_outer(defect.begin_elements()), d_end(defect.end_elements());
                        typename DenseVector<float>::ElementIterator d_inner(inner_defect.begin_elements());
                        while(d_outer != d_end )
                        {
                            *d_inner = float(*d_outer);
                            ++d_inner; ++d_outer;
                        }

                        DenseVector<float> defect_result(inner_defect.copy());
                        ConjugateGradients<Tag_, methods::JAC>::value(inner_system, inner_defect, defect_result, eps_inner);
                        inner_defect = defect_result;

                        //reconvert
                        convert(defect, inner_defect);

                        ///Update solution:
                        Scale<Tag_>::value(defect, alpha);
                        Sum<Tag_>::value(x_actual, defect);

                        defect = Product<Tag_>::value(system_matrix, x_actual);
                        Difference<Tag_>::value(defect, right_hand_side);

                        alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                        ++iter_number;



                    }
                    while(alpha < eps_outer*initial_defectnorm);

                    return x_actual;
                }

           ///sparse system:
            template<typename DT1_>
            static DenseVector<DT1_> value(SparseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                CONTEXT("When solving sparse LES with iterative refinement (PCG::JAC):");

                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                SparseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size());
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                Difference<Tag_>::value(defect, right_hand_side);

                Scale<Tag_>::value(defect, DT1_(-1.));
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Do conversion of system matrix once
                typename SparseMatrix<DT1_>::NonZeroElementIterator i_outer(system_matrix.begin_non_zero_elements()), i_end(system_matrix.end_non_zero_elements());
                typename SparseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                while(i_outer != i_end)
                {
                    while(i_inner.index() != i_outer.index())
                    {
                        ++i_inner;
                    }
                    *i_inner = float(*i_outer);
                    ++i_outer;
                }

                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(defect, DT1_(1./alpha));
                        }
                        else
                        {
                            Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }
                    DenseVector<float> defect_result(inner_defect.copy());
                    ConjugateGradients<Tag_, methods::JAC>::value(inner_system, inner_defect, defect_result, eps_inner);
                    inner_defect = defect_result;

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = DT1_(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    Scale<Tag_>::value(defect, alpha);
                    Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }


    };

    /**
     * \brief Solution of LES with Iterative Refinement using Jacobi as inner solver.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
        struct IterativeRefinement<methods::JAC, Tag_>
        {
            public:
                template<typename DT1_>
                    static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
                    {

                        CONTEXT("When solving dense LES with iterative refinement (Jacobi):");
                        DT1_ alpha = DT1_(1.0);
                        DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                        //inner allocation:
                        DenseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size(), float(0));
                        DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                        ///Initialize defect and its norm
                        //TODO: use libla/residual(?)

                        DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                        Difference<Tag_>::value(defect, right_hand_side);
                        Scale<Tag_>::value(defect, DT1_(-1.));

                        DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                        unsigned long iter_number = 0;
                        typename DenseMatrix<DT1_>::ConstElementIterator i_outer(system_matrix.begin_elements()), i_end(system_matrix.end_elements());
                        typename DenseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                        while(i_outer != i_end)
                        {
                            *i_inner = float(*i_outer);
                            ++i_inner; ++i_outer;
                        }

                        ///Main loop:
                        do
                        {
                            ///Scale defect:
                            if(iter_number != 0)
                            {
                                if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                                {
                                    Scale<Tag_>::value(defect, DT1_(1./alpha));
                                }
                                else
                                {
                                    Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                                }
                            }

                            ///Do conversion and solve inner system:
                            typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                            typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                            while(j_outer != j_end )
                            {
                                *j_inner = float(*j_outer);
                                ++j_inner; ++j_outer;
                            }

                            inner_defect = Jacobi<Tag_>::value(inner_system, inner_defect, eps_inner);

                            typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                            typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                            while(b_outer != b_end )
                            {
                                *b_outer = DT1_(*b_inner);
                                ++b_inner; ++b_outer;
                            }

                            ///Update solution:
                            Scale<Tag_>::value(defect, alpha);
                            Sum<Tag_>::value(x_actual, defect);

                            defect = Product<Tag_>::value(system_matrix, x_actual);
                            Difference<Tag_>::value(defect, right_hand_side);

                            alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                            ++iter_number;
                        }
                        while(alpha < eps_outer*initial_defectnorm);

                        return x_actual;
                    }
            ///Banded System:
            template<typename DT1_>
                static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
                {
                    CONTEXT("When solving banded LES with iterative refinement (CG):");

                    DT1_ alpha = DT1_(1.0);
                    DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                    //inner allocation:
                    BandedMatrix<float> inner_system(right_hand_side.size());
                    DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                    ///Initialize defect and its norm
                    //TODO: use libla/residual(?)

                    DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    Scale<Tag_>::value(defect, DT1_(-1.));
                    DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                    unsigned long iter_number = 0;
                    //Do conversion of system matrix once:
                    convert(inner_system, system_matrix);

                    ///Main loop:
                    do
                    {
                        ///Scale defect:
                        if(iter_number != 0)
                        {
                            if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                            {
                                Scale<Tag_>::value(defect, DT1_(1./alpha));
                            }
                            else
                            {
                                Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                            }
                        }

                        ///Do conversion and solve inner system:
                        typename DenseVector<DT1_>::ConstElementIterator d_outer(defect.begin_elements()), d_end(defect.end_elements());
                        typename DenseVector<float>::ElementIterator d_inner(inner_defect.begin_elements());
                        while(d_outer != d_end )
                        {
                            *d_inner = float(*d_outer);
                            ++d_inner; ++d_outer;
                        }

                        inner_defect = Jacobi<Tag_>::value(inner_system, inner_defect, eps_inner);
                        //reconvert
                        convert(defect, inner_defect);

                        ///Update solution:
                        Scale<Tag_>::value(defect, alpha);
                        Sum<Tag_>::value(x_actual, defect);

                        defect = Product<Tag_>::value(system_matrix, x_actual);
                        Difference<Tag_>::value(defect, right_hand_side);

                        alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                        ++iter_number;



                    }
                    while(alpha < eps_outer*initial_defectnorm);

                    return x_actual;
                }
            ///sparse system:
            template<typename DT1_>
            static DenseVector<DT1_> value(SparseMatrix<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, double eps_outer, double eps_inner)
            {
                CONTEXT("When solving sparse LES with iterative refinement (JAC):");

                DT1_ alpha = DT1_(1.0);
                DenseVector<DT1_> x_actual(right_hand_side.size(), DT1_(0));
                //inner allocation:
                SparseMatrix<float> inner_system(right_hand_side.size(), right_hand_side.size());
                DenseVector<float> inner_defect(right_hand_side.size(), float(0));

                ///Initialize defect and its norm
                //TODO: use libla/residual(?)

                DenseVector<DT1_> defect = Product<Tag_>::value(system_matrix, x_actual);
                Difference<Tag_>::value(defect, right_hand_side);

                Scale<Tag_>::value(defect, DT1_(-1.));
                DT1_ initial_defectnorm = Norm<vnt_l_two, false, Tag_>::value(defect);

                unsigned long iter_number = 0;

                ///Do conversion of system matrix once
                typename SparseMatrix<DT1_>::NonZeroElementIterator i_outer(system_matrix.begin_non_zero_elements()), i_end(system_matrix.end_non_zero_elements());
                typename SparseMatrix<float>::ElementIterator i_inner(inner_system.begin_elements());
                while(i_outer != i_end)
                {
                    while(i_inner.index() != i_outer.index())
                    {
                        ++i_inner;
                    }
                    *i_inner = float(*i_outer);
                    ++i_outer;
                }

                ///Main loop:
                do
                {
                    ///Scale defect:
                    if(iter_number != 0)
                    {
                        if(fabs(alpha) > std::numeric_limits<double>::epsilon())
                        {
                            Scale<Tag_>::value(defect, DT1_(1./alpha));
                        }
                        else
                        {
                            Scale<Tag_>::value(defect, DT1_(1./ std::numeric_limits<double>::epsilon()));
                        }
                    }

                    ///Do conversion and solve inner system:
                    typename DenseVector<DT1_>::ConstElementIterator j_outer(defect.begin_elements()), j_end(defect.end_elements());
                    typename DenseVector<float>::ElementIterator j_inner(inner_defect.begin_elements());
                    while(j_outer != j_end )
                    {
                        *j_inner = float(*j_outer);
                        ++j_inner; ++j_outer;
                    }

                    inner_defect = Jacobi<Tag_>::value(inner_system, inner_defect, eps_inner);

                    typename DenseVector<DT1_>::ElementIterator b_outer(defect.begin_elements()), b_end(defect.end_elements());
                    typename DenseVector<float>::ConstElementIterator b_inner(inner_defect.begin_elements());
                    while(b_outer != b_end )
                    {
                        *b_outer = DT1_(*b_inner);
                        ++b_inner; ++b_outer;
                    }

                    ///Update solution:
                    Scale<Tag_>::value(defect, alpha);
                    Sum<Tag_>::value(x_actual, defect);

                    defect = Product<Tag_>::value(system_matrix, x_actual);
                    Difference<Tag_>::value(defect, right_hand_side);

                    alpha = Norm<vnt_l_two, false, Tag_>::value(defect);
                    ++iter_number;
                }
                while(alpha < eps_outer*initial_defectnorm);

                return x_actual;
            }


        };

}

#endif
