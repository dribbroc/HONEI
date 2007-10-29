/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef LIBMATH_GUARD_JACOBI_HH
#define LIBMATH_GUARD_JACOBI_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/product.hh>
#include <libla/sum.hh>
#include <libla/difference.hh>
#include <libla/dot_product.hh>
#include <libla/scale.hh>
#include <libla/norm.hh>
#include <iostream>

/**
 * \file
 *
 * Templatized definition and implementation of a Jacobi Solver.
 *
 * \ingroup grpoperations
 */

using namespace std;

namespace honei
{
    /**
     * \brief Solution of LES with Jacobi method.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    template <typename Tag_=tags::CPU>
    struct Jacobi
    {
        private:
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag, DenseVector<DT1_> & diag_inverted, DenseMatrix<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<Tag_>::value(difference, former_result.copy());

                DenseVector<DT1_> temp2(right_hand_side.copy());

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2.copy();
            }
        public:
            /**
            * \brief Returns solution of LES with the Jacobi method given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param iter_number The fixed number of iterations.
            *
            */

            /// \{
            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side,long iter_number)
            {
                CONTEXT("When solving dense linear system with Jacobi (fixed # iterations):");
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                DenseMatrix<DT1_> difference(*system_matrix.copy());
                ///Create Diagonal, invert, compute difference on the fly.
                for(unsigned long i =0; i < diag.size(); ++i)
                {

                    diag[i] = system_matrix[i][i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                    difference[i][i] = DT1_(0);
                }

                DenseVector<DT1_> x(right_hand_side.copy());

                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                }
                return x;

            }

    };
}
#endif
