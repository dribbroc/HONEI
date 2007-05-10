#include <iostream>
#include <math.h>
#include "vector.hh"
#include "scalar_product.hh"

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
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

namespace pg512
{
 /** Norm is the base-class for all supported types of norms,
     * in the standard unextended case l1- and l2-norm as well as maximum-norm
     **/

 /**
     * A vector_norm_type is a template tag-parameter for the Norm template.
     * It governs the type of the norm that shall be computed.
     **/

    enum vector_norm_type {
        vnt_max = 0,
        vnt_l_one = 1,
        vnt_l_two = 2
    /// Extend here if higher L^k-norms needed
    };

    /// generic norm template - as default l2-norm is assumed
    template <typename DataType_, vector_norm_type NormType_ = vnt_l_two , bool root = false> struct VectorNorm
    {
        static DataType_ value(const Vector<DataType_> & vector)
        {
            DataType_ result(0);
            result = VectorNorm<DataType_, NormType_, root>::value();
            return result;
        }
    };

    /// partial specialisation: maximum-norm (root is ignored)
     template <typename DataType_, bool root> struct VectorNorm<DataType_, vnt_max, root>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
         for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
                {
                    if (abs(vector[l]) > result)
                    {
                        result = abs(vector[l]);
                    }
                }
                return result;
        }
    };


    /// partial specialisation: L1-norm (root is ignored)
    template <typename DataType_, bool root> struct VectorNorm<DataType_, vnt_l_one, root>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
                {
                    result += abs(vector[l]);
                }
            return result;
        }
    };

    /// total specialisation: L2-norm ^ 2
    template <typename DataType_> struct VectorNorm<DataType_, vnt_l_two, false>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
            result = VectorNorm<DataType_, vnt_l_two, false>::value();
            return result;
        }
    };

    /// total specialisation: normal L2-norm (with sqrt)
    template <typename DataType_> struct VectorNorm<DataType_, vnt_l_two, true>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
            result = VectorNorm<DataType_, vnt_l_two, true>::value();
            return result;
        }
    };

/// partial specialisation: norm with sqrt (allows L^k-norm)
    template <typename DataType_, vector_norm_type NormType_> struct VectorNorm<DataType_, NormType_, true>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result = VectorNorm<DataType_, NormType_, false>::value();
            for (int i=0; i < NormType_; ++i)
            /// if maximum-norm chosen, for is not executed because of NormType_ = 0 | assert NormType_ > 1 ?
            {
                result = sqrt(result);
            }
            return result;
        }
    };

/* Unnessacary specialisations... !?


 /// partial specialisation: norm without sqrt
    template <typename DataType_, vector_norm_type NormType_> struct VectorNorm<DataType_, NormType_, false>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
            result = VectorNorm<DataType_, NormType_, false>::value();
            return result;
        }
    };


    /// partial specialisation: L2-norm
       template <typename DataType_, bool root> struct VectorNorm<DataType_, vnt_l_two, root>
    {
        static DataType_ value (const Vector<DataType_> & vector) {
            DataType_ result(0);
            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
                {
                    result = ScalarProduct<DataType_>::value(vector, vector);
                }
            return result;
        }
    };
*/


}
