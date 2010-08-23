/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
#ifndef LIBMATH_GUARD_SQRT_HH
#define LIBMATH_GUARD_SQRT_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/util/tags.hh>

#include <cmath>

namespace honei
{
    /**
     * \brief Square root computation of the elements of the given entity.
     *
     * Sqrt is the template for the square root of the elements
     * \f[
     *     \texttt{Sqrt}(a): \quad a[i] \leftarrow \sqrt{a[i]},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grpmathoperations
     */
    template <typename Tag_ = tags::CPU> struct Sqrt
    {
        /**
         * \brief Returns the square root values of all of an entity's elements.
         *
         * \param x The entity whose elements' square root values shall be computed.
         *
         * \retval x Will modify the entity x and return it.
         */
        template <typename DT_> static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & x)
        {
            for (typename DenseVectorContinuousBase<DT_>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = std::sqrt(*i);
            }

            return x;
        }
    };

    /**
     * \brief Square root computation of the elements of the given entity.
     *
     * Sqrt is the template for the square root of the elements
     * \f[
     *     \texttt{Sqrt}(a): \quad a[i] \leftarrow \sqrt{a[i]},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grpmathoperations
     */
    template <> struct Sqrt<tags::CPU::SSE>
    {
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x);
        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x);
    };

    /**
     * \brief Square root computation of the elements of the given entity.
     *
     * Sqrt is the template for the square root of the elements
     * \f[
     *     \texttt{Sqrt}(a): \quad a[i] \leftarrow \sqrt{a[i]},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grpmathoperations
     */
    template <> struct Sqrt<tags::Cell>
    {
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x);
    };
}

#endif
