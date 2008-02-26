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
#define LIBSWE_GUARD_VOLUME_HH 1

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

#include <list>
#include <iostream>
using namespace std;
namespace honei
{
    class Volume
    {
        protected:
            double _size;
        public:

            double size()
            {
                return this->_size;
            }

            virtual void value(void)
            {
            }
    };

    class VolumeList
    {
        private:
            list<Volume *> object_list;

        public:

            virtual void convex_hull()
            {
                for(list<Volume *>::iterator i(object_list.begin()); i != object_list.end(); ++i)
                {
                    Volume * current = (*i);
                    current->value();
                }
            }

            void insert(Volume * obj)
            {
                this->object_list.push_back(obj);
            }
    };

    template<typename DataType_>
    class Cylinder : public Volume
    {
        public:
            virtual void value()
            {
                this->value(*_height, _h, _grid_x, _grid_y);
            }

            Cylinder(DenseMatrix<DataType_> & height, DataType_ h, signed long g_w, signed long g_h)
            {
                this->_height = &height;
                this->_h = h;
                this->_grid_x = g_w;
                this->_grid_y = g_h;
            }

        private:
            DenseMatrix<DataType_> * _height;
            DataType_ _h;
            signed long _grid_x;
            signed long _grid_y;

            void value(DenseMatrix<DataType_> & height, DataType_ h, signed long grid_x, signed long grid_y)
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

    template<typename DataType_>
    class Cuboid : public Volume
    {
        public:
            virtual void value()
            {
                this->value(*_height, _l, _w, _h, _grid_x, _grid_y);
            }

            Cuboid(DenseMatrix<DataType_> & height, DataType_ l, DataType_ w, DataType_ h, signed long g_w, signed long g_h)
            {
                this->_height = &height;
                this->_h = h;
                this->_w = w;
                this->_l = l;
                this->_grid_x = g_w;
                this->_grid_y = g_h;
                //normalized volume:
                this->_size = _w * _h * _l;
            }

        private:
            DenseMatrix<DataType_> * _height;
            DataType_ _h;
            DataType_ _w;
            DataType_ _l;
            signed long _grid_x;
            signed long _grid_y;

            void value(DenseMatrix<DataType_> & height, DataType_ l, DataType_ w, DataType_ h, signed long grid_x, signed long grid_y)
            {
                for(signed long i(grid_y); i < grid_y + l; ++i)
                {
                    if(i >= 0 && i < height.rows())
                    {
                        for(signed long j(grid_x); j < grid_x + w; ++j)
                        {
                            if(j >= 0 && j < height.columns())
                            {
                                height[i][j] = h;
                            }
                        }
                    }
                }
            }
    };

}
#endif
