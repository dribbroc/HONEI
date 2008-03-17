/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/libla/dense_matrix.hh>
#include <honei/libutil/tags.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Trace;

    /**
     * \brief Returns the trace of the given matrix.
     *
     * Trace is the template for the computation of the matrix's trace,
     * \f[
     *     \texttt{Trace}(a): \quad r \leftarrow \vert a(i, i) \vert.
     * \f]
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     */
    template <> struct Trace<tags::CPU>
    {
        /**
         * \{
         *
         * \brief Returns the trace of a given matrix.
         *
         * \param a The matrix whose trace shall be computed.
         *
         * \retval x Will return a scalar.
         */

        template <typename DT_>
        static DT_ value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When calculating the trace of a DenseMatrix:");

            DT_ result(0.0);

            for (unsigned long i(0), i_end(std::min(a.columns(), a.rows())) ; i < i_end ; ++i)
            {
                result += a(i, i);
            }

            return result;
        }

        /// \}
    };
}

#endif
