/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#ifndef LIBMATH_GUARD_INTERPOLATION_HH
#define LIBMATH_GUARD_INTERPOLATION_HH 1

#include <libla/dense_vector.hh>
#include <libla/difference.hh>
#include <libla/norm.hh>
#include <algorithm>

using namespace std;

namespace honei
{
    namespace interpolation_methods
    {
        class LINEAR
        {
        };
    }
    /**
     * \brief Interpolation in space with arbitrary dimensions.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     */
    template<typename Tag_, typename Method_>
    struct Interpolation
    {
    };

    template <typename Tag_>
    struct Interpolation<Tag_, interpolation_methods::LINEAR>
    {
        public:
            /**
             * \brief Interpolation in space: 2D (bilinear).
             * \param delta_x The stepsize in x direction.
             * \param delta_y The stepsize in y direction.
             * \param height The scalarfield representing the function.
             * \param x The x value in parameter space.
             * \param y The y value in parameter space.
             */
            template<typename ResPrec_>
            static inline ResPrec_ value(const ResPrec_ delta_x, const ResPrec_ delta_y, const DenseMatrix<ResPrec_>& height, const ResPrec_ x, const ResPrec_ y)
            {
                ///Compute the lower left vertex of the cell, in which the vector (x y)^T falls in
                ///parameter space:

                ResPrec_ x_t(0);
                ResPrec_ y_t(0);
                while(x_t + delta_x < x)
                {
                    x_t += delta_x;
                }
                while(y_t + delta_y < y)
                {
                    y_t += delta_y;
                }

                unsigned long i(x_t/delta_x);
                unsigned long j(y_t/delta_y);

                ///Perform bilinear interpolation:
                /*ResPrec_ result = height[i][j] * ( ResPrec_(1) - (y - y_t ))* (ResPrec_(1) - (x -x_t)) +
                                  height[i+1][j] * (y - y_t ) * (ResPrec_(1) - (y - y_t) ) +
                                  height[i][j+1] * (x -x_t) * (ResPrec_(1) - (y - y_t) ) +
                                  height[i+1][j+1] * (y - y_t ) * (x -x_t);*/

                ResPrec_ l_1 = (x - x_t) * (height[i][j+1] - height[i][j]) + height[i][j];
                ResPrec_ l_2 = (x - x_t) * (height[i+1][j+1] - height[i+1][j]) + height[i+1][j];

                ResPrec_ result = (y - y_t) * (l_2 - l_1) + l_1;

                return result;
            }
    };
}

#endif
