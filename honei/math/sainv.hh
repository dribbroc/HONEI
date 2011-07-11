/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#ifndef LIBMATH_GUARD_SAINV_HH
#define LIBMATH_GUARD_SAINV_HH 1

#include <honei/util/tags.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/product.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/difference.hh>
#include <honei/math/transposition.hh>


#include <iostream>
#include <vector>


namespace honei
{
    template<typename Tag_> struct SAINV
    {
        template <typename DT_>
        static SparseMatrix<DT_> value(const SparseMatrix<DT_> & A)
        {
            // z holds the row vector z_i
            SparseMatrix<DT_> z(A.rows(), A.columns());
            DenseVector<DT_> p(z.rows(), DT_(0));

            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                z(i, i, DT_(1));
            }

            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                SparseVector<DT_> v = Product<Tag_>::value(A, z[i]);

                for(unsigned long j(i) ; j < z.rows() ; ++j)
                {
                    p[j] = DotProduct<Tag_>::value(v, z[j]);
                }

                if (i == z.rows() - 1)
                    break;

                for(unsigned long j(i + 1) ; j < z.rows() ; ++j)
                {
                    DT_ alpha = p[j] / p[i];
                    SparseVector<DT_> z_tmp(z[i].copy());
                    Scale<Tag_>::value(z_tmp, alpha);
                    //Difference<Tag_>::value(z[j], z_tmp);
                    for (unsigned long x(0) ; x < z_tmp.size() ; ++x)
                    {
                        DT_ diff = ((const SparseMatrix<DT_>)z)[j][x] - ((const SparseVector<DT_>)z_tmp)[x];
                        if (fabs(diff) > 1e-1)
                            z[j][x] = diff;
                    }
                }
            }

            // z is the lower triangular matrix(Z^T); z_t is the upper triangular matrix (Z), z_d is the inverted diagonal matrix

            SparseMatrix<DT_> z_t(z.rows(), z.columns());
            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                z_t(i, i, DT_(1));
            }

            for (unsigned long i(0) ; i < z.rows() ; ++i)
                for (unsigned long j(0) ; j < i ; ++j)
                    if (((const SparseMatrix<DT_>)z)[i][j] != DT_(0))
                            z_t[j][i] = z[i][j];



            SparseMatrix<DT_> z_d(z.rows(), z.columns());
            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                z_d(i, i, DT_(1)/p[i]);
            }

            SparseMatrix<DT_> z_temp = Product<Tag_>::value(z_d, z);
            SparseMatrix<DT_> z_r = Product<Tag_>::value(z_t, z_temp);
            return z_r;
        }
    };
}
#endif
