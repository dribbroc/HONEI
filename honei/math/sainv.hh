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
#include <honei/la/sparse_matrix_ell.hh>

#include <iostream>
#include <vector>
#include <algorithm>

namespace
{
    template <typename DT_>
    struct matrix_tupple
    {
        unsigned long row;
        unsigned long column;
        DT_ value;

        static bool compare (matrix_tupple<DT_> x, matrix_tupple<DT_> y)
        {
            return (x.value) > (y.value);
        }
    };
}

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
            SparseMatrixELL<DT_> Aell(A);

            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                z(i, i, DT_(1));
            }

            DenseVector<DT_> v(Aell.columns());
            for(unsigned long i(0) ; i < z.rows() ; ++i)
            {
                DenseVector<DT_> zi(((const SparseMatrix<DT_>)z)[i]);
                Product<Tag_>::value(v, Aell, zi);
                /*SparseVector<DT_> v2(v.size(), 1);
                for (unsigned long vi(0) ; vi < v.size() ; ++vi)
                    if (v[vi] != DT_(0))
                        v2[vi] = v[vi];*/

                /*for(unsigned long j(i) ; j < z.rows() ; ++j)
                {
                    p[j] = DotProduct<tags::CPU>::value(v2, z[j]);
                }*/
                const DT_ * vele(v.elements());
                DT_ * pele(p.elements() + i);
                for(unsigned long j(i) ; j < z.rows() ; ++j, ++pele)
                {
                    const DT_ * ele(z[j].elements());
                    const unsigned long * ind(z[j].indices());
                    *pele = DT_(0);
                    for (unsigned long k(0) ; k < (z[j]).used_elements() ; ++k)
                    {
                        *pele += ele[k] * vele[ind[k]];
                    }

                }

                if (i == z.rows() - 1)
                    break;

                for(unsigned long j(i + 1) ; j < z.rows() ; ++j)
                {
                    //DT_ alpha = p[j] / p[i];
                    DT_ alpha = (fabs(p[i]) > std::numeric_limits<DT_>::epsilon()) ? p[j] / p[i] : p[j] / std::numeric_limits<DT_>::epsilon();
                    if (fabs(alpha) > tolerance)
                    {
                        //TODO integrate dropping strategy in difference and insert (virtual) zeros if dropping occurs
                        typename SparseVector<DT_>::NonZeroElementIterator l(z[j].begin_non_zero_elements());
                        typename SparseVector<DT_>::NonZeroConstElementIterator r(z[i].begin_non_zero_elements());
                        typename SparseVector<DT_>::NonZeroConstElementIterator r_end(z[i].end_non_zero_elements());
                        for ( ; r != r_end ; )
                        {
                            if (r.index() < l.index())
                            {
                                if(fabs(*r * alpha) > tolerance)
                                {
                                    z[j][r.index()] = -(*r) * alpha;
                                }
                                ++r;
                            }
                            else if (l.index() < r.index())
                            {
                                ++l;
                            }
                            else
                            {
                                //TODO BUG!: if dropping occurs: do not insert zero but delete the entry at z[j][x] if any exists
                                *l = (*l - (*r * alpha)) * DT_(fabs(*l - (*r * alpha)) > tolerance);
                                ++l; ++r;
                            }
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

            // TODO Tag_ nicht hardverdrahten
            SparseMatrix<DT_> z_temp = Product<tags::CPU>::value(z_d, z);
            SparseMatrix<DT_> z_r = Product<tags::CPU>::value(z_t, z_temp);

            //return z_r;

            //POST FILTERING
            SparseMatrix<DT_> result(A.rows(), A.columns());
            for (unsigned long i(0) ; i < A.rows() ; ++i)
            {
                std::vector<matrix_tupple<DT_> > elements;

                typename SparseVector<DT_>::NonZeroConstElementIterator r(z_r[i].begin_non_zero_elements());
                typename SparseVector<DT_>::NonZeroConstElementIterator r_end(z_r[i].end_non_zero_elements());
                for ( ; r != r_end ; ++r)
                {
                    matrix_tupple<DT_> temp;
                    temp.row = i;
                    temp.column = r.index();
                    temp.value = *r;
                    elements.push_back(temp);
                }

                unsigned long target(std::max(A.used_elements() / A.rows(), A[i].used_elements()));
                if (elements.size() > target)
                {
                    // sort from higher to lower
                    std::sort(elements.begin(), elements.end(), matrix_tupple<DT_>::compare);
                    // remove smalles values
                    elements.resize(target);
                }

                for (unsigned long i(0) ; i < elements.size() ; ++i)
                {
                    result(elements.at(i).row, elements.at(i).column, elements.at(i).value);
                }
            }

            return result;
        }
    };
}
#endif
