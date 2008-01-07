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

        class NN
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
            static inline ResPrec_ value(const ResPrec_ delta_x, const ResPrec_ delta_y, const DenseMatrix<ResPrec_>& height, ResPrec_ x, ResPrec_ y)
            {
                ///Compute the lower left vertex of the cell, in which the vector (x y)^T falls in
                ///parameter space:
                ResPrec_ x_t_1(x / delta_x);
                ResPrec_ y_t_1(y / delta_y);

                ///Catch GAMMA case (negative):
                if(x_t_1 < ResPrec_(0))
                {
                    x_t_1 = ResPrec_(0);
                    x = ResPrec_(0);
                }
                if(y_t_1 < ResPrec_(0))
                {
                    y_t_1 = ResPrec_(0);
                    y = ResPrec_(0);
                }

                unsigned long x_t = (unsigned long)x_t_1;
                unsigned long y_t = (unsigned long)y_t_1;
                unsigned long j = (x_t);
                unsigned long i = (y_t);
                ResPrec_ l_1, l_2;

                ///Perform bilinear interpolation:
                if(i < height.rows() - 1 && j < height.columns() - 1)
                {
                    l_1 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[i][j+1] - height[i][j]) + height[i][j];
                    l_2 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[i+1][j+1] - height[i+1][j]) + height[i+1][j];
                }
                else if(i >= height.rows() - 1 && j >= height.columns() - 1)
                {
                    /*l_1 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[height.rows() - 1][height.columns() - 1] - height[height.rows() - 1][height.columns() - 1]) + height[height.rows() - 1][height.columns() - 1];
                    l_2 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[height.rows() - 1][height.columns() - 1] - height[height.rows() - 1][height.columns() - 1]) + height[height.rows() - 1][height.columns() - 1];*/
                    if(i > j)
                    {
                        return height[i][height.columns() - 1];
                    }
                    else
                    {
                        return height[height.rows() - 1][j];
                    }
                }
                else if(i >= height.rows() - 1)
                {
                    /*l_1 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[height.rows() - 1][j+1] - height[height.rows() - 1][j]) + height[height.rows() - 1][j];
                    l_2 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[height.rows() - 1][j+1] - height[height.rows() - 1][j]) + height[height.rows() - 1][j];*/
                    return height[height.rows()][j];
                }
                else if(j >= height.columns() - 1)
                {
                    /*l_1 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[i][height.columns() - 1] - height[i][height.columns() - 1]) + height[i][height.columns() - 1];
                    l_2 = (x - ResPrec_(j * delta_x))/(ResPrec_(((j + 1) * delta_x) - (j * delta_x))) * (height[i + 1][height.columns() - 1] - height[i + 1][height.columns() - 1]) + height[i+1][height.columns() - 1];*/
                    return height[i][height.columns() - 1];

                }
                ResPrec_ result = (y - ResPrec_(i * delta_y))/(ResPrec_(((i + 1) * delta_y) - (i * delta_y))) * (l_2 - l_1) + l_1;
                return result;
            }
    };

    template <typename Tag_>
        struct Interpolation<Tag_, interpolation_methods::NN>
        {
            public:
                /**
                 * \brief Interpolation in space: 2D (NEAREST NEIGHBOUR).
                 * \param delta_x The stepsize in x direction.
                 * \param delta_y The stepsize in y direction.
                 * \param height The scalarfield representing the function.
                 * \param x The x value in parameter space.
                 * \param y The y value in parameter space.
                 */
                template<typename ResPrec_>
                static inline ResPrec_ value(const ResPrec_ delta_x, const ResPrec_ delta_y, const DenseMatrix<ResPrec_>& height, ResPrec_ x, ResPrec_ y)
                {
                    ///Compute the lower left vertex of the cell, in which the vector (x y)^T falls in
                    ///parameter space:
                    ResPrec_ x_t_1(x / delta_x);
                    ResPrec_ y_t_1(y / delta_y);
                    unsigned long x_t = (unsigned long)x_t_1;
                    unsigned long y_t = (unsigned long)y_t_1;
                    unsigned long j = (x_t);
                    unsigned long i = (y_t);
                    ///Catch GAMMA case (negative):
                    if(x_t_1 < ResPrec_(0))
                    {
                        x_t_1 = ResPrec_(0);
                        x = ResPrec_(0);
                    }
                    if(y_t_1 < ResPrec_(0))
                    {
                        y_t_1 = ResPrec_(0);
                        y = ResPrec_(0);
                    }

                    unsigned long nearest_x;
                    unsigned long nearest_y;
                    if(x < ResPrec_(j * delta_x) + ResPrec_(delta_x/2))
                    {
                        nearest_x = j;
                    }
                    else
                    {
                        nearest_x = j + 1;
                    }

                    if(y < ResPrec_(i * delta_y) + ResPrec_(delta_y/2))
                    {
                        nearest_y = i;
                    }
                    else
                    {
                        nearest_y = i + 1;
                    }

                    if(nearest_x < height.columns() - 1 && nearest_y < height.rows() - 1)
                    {
                        return height[nearest_y][nearest_x];
                    }
                    else if(nearest_x >= height.columns() - 1 && nearest_y >= height.rows() - 1)
                    {
                        return height[height.rows() - 1][height.columns() -1];
                    }
                    else if (nearest_x >= height.columns() - 1)
                    {
                        return height[nearest_y][height.columns() - 1];
                    }
                    else if (nearest_y >= height.rows() - 1)
                    {
                        return height[height.rows() - 1][nearest_x];
                    }
                }
        };
}

#endif
