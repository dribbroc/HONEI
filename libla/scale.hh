/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_SCALE_HH
#define LIBLA_GUARD_SCALE_HH 1

#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libutil/tags.hh>

/**
 * \file
 *
 * Templatized definitions of operation Scale.
 *
 * \ingroup grpoperations
 */
namespace honei
{
    /**
     * \brief Result of scaling an entity by a scalar factor.
     *
     * Scale is the class template for the operation
     * \f[
     *     Scale(a, x): \quad x \leftarrow a * x,
     * \f]
     * which yields the former entity scaled by a scalar factor.
     *
     * \ingroup grpoperations
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_ = tags::CPU> struct Scale
    {
        /**
         * \brief Returns the product of a scalar factor and a given entity.
         *
         * \param a The scalar factor.
         * \param x The entity that shall be scaled.
         *
         * \retval x Will modify the entity x and return it.
         */

        /// \{

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const DT1_ a, DenseMatrix<DT2_> & x)
        {
            for (typename MutableMatrix<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(const DT1_ a, SparseMatrix<DT2_> & x)
        {
            for (typename MutableMatrix<DT2_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(const DT1_ a, BandedMatrix<DT2_> & x)
        {
            for (typename BandedMatrix<DT2_>::VectorIterator l(x.begin_bands()),
                    l_end(x.end_bands()) ; l != l_end ; ++l)
            {
                *l = Scale<>::value(a, *l);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT2_> & value(const DT1_ a, DenseVector<DT2_> & x)
        {
            for (typename Vector<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> & value(const DT1_ a, SparseVector<DT2_> & x)
        {
            for (typename Vector<DT2_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        /// \}
    };
}
#endif
