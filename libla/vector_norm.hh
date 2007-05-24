/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_VECTOR_NORM_HH
#define LIBLA_GUARD_VECTOR_NORM_HH 1

#include <libla/tags.hh>
#include <libla/vector.hh>
#include <libla/scalar_product.hh>

#include <cmath>
#include <iostream>

/**
 * \file
 *
 * Templatized definitions of vector norms.<br/>
 *
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * A VectorNormType is a template tag parameter for the VectorNorm class
     * template. It governs the type of the (mathematical) norm that shall be
     * computed.
     *
     * \ingroup grpvectoroperations
     **/
    enum VectorNormType {
        vnt_max = 0,
        vnt_l_one,
        vnt_l_two,
        // Extend here if higher norms needed
    };

    /**
     * VectorNorm is the class template for all supported types of vector norms.
     * The standard norms of type L1 and L2 are supported, as well as the
     * maximum norm. As default tag the L2 vector-norm is assumed.
     *
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, VectorNormType norm_type_ = vnt_l_two,
             bool root_ = false, typename Tag_ = tags::CPU> struct VectorNorm
    {
        /**
         * Return the norm of a vector.
         *
         * \param vector Vector whose norm shall be computed.
         **/
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            DataType_ result(0);

            result = VectorNorm<DataType_, norm_type_, root_>::value(vector);

            return result;
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            DataType_ result(0);

            result = VectorNorm<DataType_, norm_type_, root_>::value(vector);

            return result;
        }
    };

    /// Partial specialisation of VectorNorm for vnt_max (maximums-norm)
    /// \ingroup grpvectoroperations
    template <typename DataType_, bool root_> struct VectorNorm<DataType_, vnt_max, root_, tags::CPU>
    {
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ;
                    l != l_end ; ++l)
            {
                if (fabs(*l) > result)
                {
                    result = fabs(*l);
                }
            }

            return result;
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_non_zero_elements()), l_end(vector.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                if (abs(*l) > result)
                {
                    result = abs(*l);
                }
            }

            return result;
        }
    };


    /// Partial specialisation of VectorNorm for vnt_l_one (l1-norm)
    /// \ingroup grpvectoroperations
    template <typename DataType_, bool root_> struct VectorNorm<DataType_, vnt_l_one, root_, tags::CPU>
    {
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ;
                    l != l_end ; ++l)
            {
                result += fabs(*l);
            }

            return result;
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_non_zero_elements()), l_end(vector.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                result += fabs(*l);
            }

            return result;
        }
    };

    /// Partial specialisation of VectorNorm for vnt_l_two (square of L2-norm)
    /// \ingroup grpvectoroperations
    template <typename DataType_> struct VectorNorm<DataType_, vnt_l_two, false, tags::CPU>
    {
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            return ScalarProduct<DataType_>::value(vector, vector);
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            return ScalarProduct<DataType_>::value(vector, vector);
        }
    };

    /// Partial specialisation of VectorNorm for vnt_l_two (L2-norm)
    /// \ingroup grpvectoroperations
    template <typename DataType_> struct VectorNorm<DataType_, vnt_l_two, true, tags::CPU>
    {
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            return sqrt(VectorNorm<DataType_, vnt_l_two, false>::value(vector));
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            return sqrt(VectorNorm<DataType_, vnt_l_two, false>::value(vector));
        }
    };

    /// Partial specialisation of VectorNorm for arbitrary norms (allows Lk-norm)
    /// \ingroup grpvectoroperations
    template <typename DataType_, VectorNormType norm_type_> struct VectorNorm<DataType_,
             norm_type_, true, tags::CPU>
    {
        static DataType_ value (const DenseVector<DataType_> & vector)
        {
            DataType_ result(0);
            unsigned int k(static_cast<unsigned int>(norm_type_));

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ;
                    l != l_end ; ++l)
            {
                if (*l != static_cast<DataType_>(0))
                {
                    result += exp(k * log(abs(*l)));
                }
            }

            return exp(log(result) / k);
        }

        static DataType_ value (const SparseVector<DataType_> & vector)
        {
            DataType_ result(0);
            unsigned int k(static_cast<unsigned int>(norm_type_));

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_non_zero_elements()), l_end(vector.end_non_zero_elements()) ;
                    l != l_end ; ++l)
            {
                result += exp(k * log(fabs(*l)));
            }

            return exp(log(result) / k);
        }
    };
}

#endif
