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

#ifndef LIBLA_GUARD_VECTOR_ABSOLUTE_HH
#define LIBLA_GUARD_VECTOR_ABSSOLUTE_HH 1

#include <libla/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>

#include <cmath>

/**
 * \file
 *
 * Templatized definitions of calculating a vector's absolute values.<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{
    /**
     * VectorAbsolute is the class template that calculates a vector's absolute
     * values.
     *
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct VectorAbsolute
    {
        /**
         * Return a vector's absolute values.
         *
         * \param vector DenseVector whose absolute values shall be computed.
         **/
        static DenseVector<DataType_> & value(DenseVector<DataType_> & vector)
        {
            for (typename Vector<DataType_>::ElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i < static_cast<DataType_>(0))
                    *i *= -1;
            }

            return vector;
        }
    };

    /// Specialisation of VectorAbsolute for float
    /// \ingroup grpvectoroperations
    template <> struct VectorAbsolute<float, tags::CPU>
    {
        /**
         * Return a vector's absolute values.
         *
         * \param vector DenseVector whose absolute values shall be computed.
         **/
        static DenseVector<float> & value(DenseVector<float> & vector);

        /**
         * Return a vector's absolute values.
         *
         * \param vector SparseVector whose absolute values shall be computed.
         **/
        static SparseVector<float> & value(SparseVector<float> & vector);
    };

    /// Specialisation of VectorAbsolute for double
    /// \ingroup grpvectoroperations
    template <> struct VectorAbsolute<double, tags::CPU>
    {
        /**
         * Return a vector's absolute values.
         *
         * \param vector DenseVector whose absolute values shall be computed.
         **/
        static DenseVector<double> & value(DenseVector<double> & vector);

        /**
         * Return a vector's absolute values.
         *
         * \param vector SparseVector whose absolute values shall be computed.
         **/
        static SparseVector<double> & value(SparseVector<double> & vector);
    };
}

#endif
