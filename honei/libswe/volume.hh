/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBSWE_GUARD_VOLUME_HH
#define LIBSWE_GUARD_VOLUME_SHOT_HH 1

/**
 * \file
 *
 * Implementation of utility classes to generate meta-volumes on the water or ground
 * surface for use with ScenarioManger.
 *
 * \ingroup grplibswe
 **/

#include <honei/libla/dense_matrix.hh>
#include <honei/libla/dense_vector.hh>

namespace volume_types
{
    class CYLINDRIC
    {
        class STENCIL;
        class BRESENHAM;
    };
}

namespace honei
{
    template<typename VolumeType_>
    struct Volume
    {
    };

    template<>
    struct Volume<volume_types::CYLINDRIC::STENCIL>
    {
        public:
            template<typename DataType_>
            static DenseMatrix<DataType_>& value(DenseMatrix<DataType_> & height, DataType_ h, signed long grid_x, signed long grid_y)
            {
                if(grid_y >= 0 && grid_y < height.rows())
                {
                    for(signed long i(grid_x - 6); i < grid_x + 6; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y][i] = h;
                        }
                    }
                }

                if(grid_y + 1 >= 0 && grid_y + 1 < height.rows())
                {
                    for(signed long i(grid_x - 6); i < grid_x + 6; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 1][i] = h;
                        }
                    }
                }

                if(grid_y - 1 >= 0 && grid_y - 1 < height.rows())
                {
                    for(signed long i(grid_x - 6); i < grid_x + 6; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 1][i] = h;
                        }
                    }
                }

                if(grid_y + 2 >= 0 && grid_y + 2 < height.rows())
                {
                    for(signed long i(grid_x - 5); i < grid_x + 5; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 2][i] = h;
                        }
                    }

                }

                if(grid_y - 2 >= 0 && grid_y - 2 < height.rows())
                {
                    for(signed long i(grid_x - 5); i < grid_x + 5; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 2][i] = h;
                        }
                    }

                }

                if(grid_y + 3 >= 0 && grid_y + 3 < height.rows())
                {
                    for(signed long i(grid_x - 5); i < grid_x + 5; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 3][i] = h;
                        }
                    }

                }

                if(grid_y - 3 >= 0 && grid_y - 3 < height.rows())
                {
                    for(signed long i(grid_x - 5); i < grid_x + 5; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 3][i] = h;
                        }
                    }

                }

                if(grid_y + 4 >= 0 && grid_y + 4 < height.rows())
                {
                    for(signed long i(grid_x - 4); i < grid_x + 4; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 4][i] = h;
                        }
                    }

                }

                if(grid_y - 4 >= 0 && grid_y - 4 < height.rows())
                {
                    for(signed long i(grid_x - 4); i < grid_x + 4; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 4][i] = h;
                        }
                    }

                }

                if(grid_y + 5 >= 0 && grid_y + 5 < height.rows())
                {
                    for(signed long i(grid_x - 3); i < grid_x + 3; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 5][i] = h;
                        }
                    }

                }

                if(grid_y - 5 >= 0 && grid_y - 5 < height.rows())
                {
                    for(signed long i(grid_x - 3); i < grid_x + 3; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 5][i] = h;
                        }
                    }

                }

                if(grid_y + 6 >= 0 && grid_y + 6 < height.rows())
                {
                    for(signed long i(grid_x - 2); i < grid_x + 2; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y + 6][i] = h;
                        }
                    }

                }

                if(grid_y - 6 >= 0 && grid_y - 6 < height.rows())
                {
                    for(signed long i(grid_x - 2); i < grid_x + 2; ++i)
                    {
                        if(i >= 0 && i < height.columns())
                        {
                            height[grid_y - 6][i] = h;
                        }
                    }

                }

            }
    };
}
#endif
