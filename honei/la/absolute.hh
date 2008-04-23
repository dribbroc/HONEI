/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_ABSOLUTE_HH
#define LIBLA_GUARD_ABSOLUTE_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/vector.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/util/tags.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Absolute;

    /**
     * \brief Returns the absolute of the given entity.
     *
     * Absolute is the template for the computation of the entity's absolute,
     * \f[
     *     \texttt{Absolute}(a): \quad a[i] \leftarrow \vert a[i] \vert.
     * \f]
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Absolute<tags::CPU>
    {
        /**
         * \name Vector absolutes
         * \{
         *
         * \brief Returns the absolute of a given vector.
         *
         * \param x The vectors whose elements' absolute values shall be computed.
         *
         * \retval x Will modify the vector x and return it.
         */

        template <typename DT_>
        static DenseVectorBase<DT_> & value(DenseVectorBase<DT_> & x)
        {
            CONTEXT("When calculating the absolute value of DenseVectorBase elements");

            for (typename Vector<DT_>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i < DT_(0))
                    *i *= DT_(-1);
            }

            return x;
        }

        template <typename DT_>
        static SparseVector<DT_> & value(SparseVector<DT_> & x);

        /// \}

        template <typename DT_>
        static inline DenseVector<DT_> & value(DenseVector<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            Absolute<>::value(temp);
            return x;
        }

        template <typename DT_>
        static inline DenseVectorSlice<DT_> & value(DenseVectorSlice<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            Absolute<>::value(temp);
            return x;
        }

        template <typename DT_>
        static inline DenseVectorRange<DT_> & value(DenseVectorRange<DT_> & x)
        {
            DenseVectorBase<DT_> & temp = x;
            Absolute<>::value(temp);
            return x;
        }
    };

    template <>
    DenseVectorBase<float> & Absolute<tags::CPU>::value(DenseVectorBase<float> & x);

    template <>
    SparseVector<float> & Absolute<tags::CPU>::value(SparseVector<float> & x);

    template <>
    DenseVectorBase<double> & Absolute<tags::CPU>::value(DenseVectorBase<double> & x);

    template <>
    SparseVector<double> & Absolute<tags::CPU>::value(SparseVector<double> & x);
}

#endif
