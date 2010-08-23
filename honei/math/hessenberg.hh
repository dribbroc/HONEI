/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@tu-dortmund.de>
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

#pragma once
#ifndef LIBMATH_GUARD_HESSENBERG_HH
#define LIBMATH_GUARD_HESSENBERG_HH 1

#include <honei/util/tags.hh>
#include <honei/la/algorithm.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_matrix_tile.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/scale.hh>
#include <honei/la/norm.hh>
#include <honei/la/product.hh>

namespace honei
{
    template<typename Tag_> struct Hessenberg
    {
    };

    /**
     * Transforms a matrix to Hessenberg form.
     *
     * The referenced containers are modified under this operation.
     *
     * \ingroup grpmatrixoperations
     */
    template <> struct Hessenberg<tags::CPU>
    {
        private:
            template <typename DT_> static inline void _step(DenseMatrix<DT_> & matrix, unsigned long index)
            {
                unsigned long size(matrix.columns() - 1 - index);

                CONTEXT("When performing " + stringify(index) + "th step of Hessenberg decomposition:");
                ASSERT(size > 0, "invalid index!");

                DenseVectorSlice<DT_> column(matrix.column(index));

                DenseVector<DT_> x(size);
                honei::copy(column.element_at(1 + index), column.end_elements(), x);

                DT_ c(Norm<vnt_l_two, true, tags::CPU>::value(x));
                if (c <= std::numeric_limits<DT_>::epsilon())
                    return;

                column[1 + index] = -c;
                honei::fill(column.element_at(2 + index), column.end_elements(), DT_(0));

                DenseVectorRange<DT_> y(matrix[index], size, 1 + index);

                DenseVector<DT_> w(x);
                w[0] += c;
                Scale<tags::CPU>::value(w, DT_(1.0) / Norm<vnt_l_two, true, tags::CPU>::value(w));

                ScaledSum<tags::CPU>::value(y, w, DT_(-2.0) * DotProduct<tags::CPU>::value(w, y));

                DenseMatrixTile<DT_> a(matrix, size, size, matrix.rows() - size, matrix.columns() - size);

                /// \todo Not yet implemented: DenseVector<DT_> q_l(Product<tags::CPU>::value(w, a))
                DenseVector<DT_> q_l(Product<tags::CPU>::value(w, a));
#if 0
                DenseVector<DT_> q_l(a[0].copy());
                Scale<tags::CPU>::value(q_l, w[0]);
                for (unsigned i(1) ; i < size ; ++i)
                {
                    ScaledSum<tags::CPU>::value(q_l, a[i], w[i]);
                }
#endif

                /// \todo Not yet implemented: DenseVector<DT_> q_r(Product<tags::CPU>::value(a, w))
                DenseVector<DT_> q_r(size);
                for (unsigned i(0) ; i < size ; ++i)
                {
                    q_r[i] = DotProduct<tags::CPU>::value(a[i], w);
                }

                DT_ r(DotProduct<tags::CPU>::value(w, q_r));

                for (unsigned i(0) ; i < size ; ++i)
                {
                    ScaledSum<tags::CPU>::value(a[i], q_l, DT_(-2.0) * w[i]);

                    ScaledSum<tags::CPU>::value(a[i], w, DT_(-2.0) * q_r[i]);

                    ScaledSum<tags::CPU>::value(a[i], w, DT_(4.0) * r * w[i]);
                }
            }

        public:
            template <typename DT_> static inline void value(DenseMatrix<DT_> & matrix)
            {
                CONTEXT("When performing Hessenberg decomposition:");

                if (matrix.rows() != matrix.columns())
                    throw MatrixIsNotSquare(matrix.rows(), matrix.columns());

                for (unsigned i(0) ; i < matrix.rows() - 2 ; ++i)
                {
                    _step(matrix, i);
                }
            }
    };
}

#endif
