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

// Based on "Robust Approximate Inverse Preconditioning for the Conjugate Gradients Method" by Benzi et al.
namespace honei
{
    template<typename Tag_> struct SAINV
    {
        template <typename DT_>
        static SparseMatrix<DT_> value(const SparseMatrix<DT_> & A, DT_ tolerance = 13e-2)
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
                    //TODO integrate dropping strategy in difference and insert (virtual) zeros if dropping occurs
                    typename SparseVector<DT_>::NonZeroElementIterator l(z[j].begin_non_zero_elements());
                    for (typename SparseVector<DT_>::NonZeroConstElementIterator r(z_tmp.begin_non_zero_elements()),
                            r_end(z_tmp.end_non_zero_elements()) ; r != r_end ; )
                    {
                        if (r.index() < l.index())
                        {
                            if(fabs(*r) > tolerance)
                                z[j][r.index()] = -(*r);
                            ++r;
                        }
                        else if (l.index() < r.index())
                        {
                            ++l;
                        }
                        else
                        {
                            //TODO BUG!: if dropping occurs: do not insert zero but delete the entry at z[j][x] if any exists
                            *l = (*l - *r) * DT_(fabs(*l - *r) > tolerance);
                            ++l; ++r;
                        }
                    }
                }
            }

            // z is the lower triangular matrix(Z^T); z_t is the upper triangular matrix (Z), z_d is the inverted diagonal matrix

            SparseMatrix<DT_> z_t(z.rows(), z.columns());
            Transposition<Tag_>::value(z, z_t);
            for(unsigned long i(0) ; i < z_t.rows() ; ++i)
            {
                z_t(i, i, DT_(1));
            }

            SparseMatrix<DT_> z_d(z.rows(), z.columns());
            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                z_d(i, i, DT_(1)/p[i]);
            }

            SparseMatrix<DT_> z_temp = Product<Tag_>::value(z_d, z);
            SparseMatrix<DT_> z_r = Product<Tag_>::value(z_t, z_temp);

            /* POST FILTERING
            unsigned long ue(0);
            for (typename SparseMatrix<DT_>::NonZeroElementIterator r(z_r.begin_non_zero_elements()),
                    r_end(z_r.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                if (fabs(*r) < tolerance / 5)
                {
                    ue++;
                    *r = DT_(0);
                }
            }
            std::cout<<ue<<std::endl;*/
            return z_r;
        }
    };
}
#endif
