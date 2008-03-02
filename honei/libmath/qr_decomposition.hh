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

#ifndef LIBMATH_GUARD_QR_DECOMPOSITION_HH
#define LIBMATH_GUARD_QR_DECOMPOSITION_HH 1

#include <honei/libutil/tags.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/scale.hh>
#include <honei/libla/scaled_sum.hh>

#include <cmath>

namespace honei
{
    template<typename Tag_> struct QRDecomposition
    {
    };

    /**
     * Transforms a matrix in Hessenberg form to upper triangular form.
     *
     * The referenced containers are modified under this operation.
     *
     * \ingroup grpmatrixoperations
     */
    template <> struct QRDecomposition<tags::CPU>
    {
        private:
            template <typename DT_> static inline void _step(DenseMatrix<DT_> & matrix, unsigned long index)
            {
                CONTEXT("When performing " + stringify(index) + "th QR decomposition step:");
                ASSERT(index < matrix.rows() - 1, "invalid index!");

                DenseVectorRange<DT_> row(matrix[index + 1]);

                if (fabs(matrix(index + 1, index)) <= std::numeric_limits<DT_>::epsilon())
                    return;

                DenseVectorRange<DT_> row_a(matrix[index]);
                DenseVector<DT_> row_a_copy(row_a.copy());

                DenseVectorRange<DT_> row_b(matrix[index + 1]);
                DenseVector<DT_> row_b_copy(row_b.copy());

                DT_ c(0.0);
                DT_ s(1.0);

                if (fabs(row_a[index]) > std::numeric_limits<DT_>::epsilon())
                {
                    DT_ t(row_b[index] / row_a[index]);
                    DT_ theta(atan(t));

                    s = sin(theta);
                    c = cos(theta);
                }

                Scale<tags::CPU>::value(row_a, c);
                ScaledSum<tags::CPU>::value(row_a, row_b_copy, s);

                Scale<tags::CPU>::value(row_b, c);
                ScaledSum<tags::CPU>::value(row_b, row_a_copy, DT_(-1.0) * s);
            }

        public:
            template <typename DT_> static inline void value(DenseMatrix<DT_> & matrix)
            {
                CONTEXT("When performing QR decomposition:");

                for (unsigned i(0) ; i < matrix.rows() - 1 ; ++i)
                {
                    _step(matrix, i);
                }
            }
    };
}

#endif
