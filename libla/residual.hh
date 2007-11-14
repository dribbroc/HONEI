 /* vim: set sw=4 sts=4 et nofoldenable : */

 /*
  * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
  *
  * This file is part of the LA C++ library. LibLa is free software;
  * you can redistribute it and/or modify it under the terms of the GNU Genera
  * Public License version 2, as published by the Free Software Foundation.
  *
  * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
  * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  * details.
  *
  * You should have received a copy of the GNU General Public License along with
  * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
  */

#ifndef LIBLA_GUARD_RESIDUAL_HH
#define LIBLA_GUARD_RESIDUAL_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/dot_product.hh>
#include <libla/matrix_error.hh>
#include <libla/sparse_matrix.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

namespace honei
{
    /**
     * \brief Residual of a system of linear equations.
     *
     * Residual is the class template for the operation
     * \f[
     *    \texttt{Residual}(b, A, x): \quad b \leftarrow b - A \cdot x,
     * \f]
     * which yields the residual or defect of a given approximative solution to
     * the system of linear equations
     * \f[
     *     A \cdot x = b.
     * \f]
     *
     * \ingroup grplaoperations
     */
    template <typename Tag_ = tags::CPU> struct Residual
    {
        /**
         * \name Residuals
         * \{
         *
         * Returns the the residual or defect of a given approximative solution
         * to a linear equation.
         *
         * \param b The vector of inhomogenous parts of the system of linear equations.
         * \param a The matrix that corresponds to the system of lineare
         *          equations.
         * \param x The vector of approximative solutions to the system of linear
         *          equations.
         *
         * \retval b Will modify b and return it.
         */

        template <typename DT1_, typename DT2_, typename VT_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & b,
                const DenseMatrix<DT2_> & a,
                const VT_ & x)
        {
            CONTEXT("When calculating residual for a DenseMatrix:");

            if (b.size() != x.size())
                throw VectorSizeDoesNotMatch(x.size(), b.size());

            if (a.rows() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.rows());

            if (! a.square())
                throw MatrixIsNotSquare(a.rows(), a.columns());

            for (typename Vector<DT1_>::ElementIterator i(b.begin_elements()),
                    i_end(b.end_elements()) ; i != i_end ; ++i)
            {
                *i -= DotProduct<>::value(a[i.index()], x);
            }

            return b;
        }

        template <typename DT1_, typename DT2_, typename VT_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & b,
                const SparseMatrix<DT2_> & a,
                const VT_ & x)
        {
            CONTEXT("When calculating residual for a SparseMatrix:");

            if (b.size() != x.size())
                throw VectorSizeDoesNotMatch(x.size(), b.size());

            if (a.rows() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.rows());

            if (! a.square())
                throw MatrixIsNotSquare(a.rows(), a.columns());

            for (typename Vector<DT1_>::ElementIterator i(b.begin_elements()),
                    i_end(b.end_elements()) ; i != i_end ; ++i)
            {
                *i -= DotProduct<>::value(a[i.index()], x);
            }

            return b;
        }

        /// \}
    };
}

#endif
