/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
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
#include <unittest/unittest.hh>
#include <iostream>
#include<honei/lbm/solid.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;

template <typename Tag_, typename DataType_>
class SolidTest :
    public TaggedTest<Tag_>
{
    public:
        SolidTest(const std::string & type) :
            TaggedTest<Tag_>("solid_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            Line<DataType_, lbm_solid_dims::D2> line_0(DataType_(1), DataType_(1), DataType_(10), DataType_(10));
            Line<DataType_, lbm_solid_dims::D2> line_1(DataType_(10), DataType_(10), DataType_(19), DataType_(1));
            Line<DataType_, lbm_solid_dims::D2> line_2(DataType_(19), DataType_(1), DataType_(1), DataType_(1));

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(3);
            tri_0.add_line(line_2);
            tri_0.add_line(line_1);
            tri_0.add_line(line_0);

            tri_0.value();

            std::cout << "Polygon with " << tri_0.get_line_count() << " lines: " << std::endl;
            for (unsigned long i(0); i < 6 ; ++i)
            {
                std::cout << "Vertex at (" << tri_0.get_x_coords()[i] << " , " << tri_0.get_y_coords()[i] << ")" << std::endl;
            }

            Polygon<DataType_, lbm_solid_dims::D2> poly_0(1000);
            for(unsigned long i(0) ; i < 1000; ++i)
            {
                Line<DataType_, lbm_solid_dims::D2> line_i(DataType_(i), DataType_(i), DataType_(i+1), DataType_(i+1));
                poly_0.add_line(line_i);
            }

            poly_0.value();
            std::cout << "Polygon with " << poly_0.get_line_count() << " lines: " << std::endl;
            for (unsigned long i(0); i < 2000 ; ++i)
            {
                std::cout << "Vertex at (" << poly_0.get_x_coords()[i] << " , " << poly_0.get_y_coords()[i] << ")" << std::endl;
            }
            //Rasterization test:
            DenseMatrix<bool> target(20, 20, false);
            ScanConversion<Tag_>::value(tri_0, target, DataType_(1), DataType_(1));

            std::cout << target;
        }
};
SolidTest<tags::CPU, float> solidtest_float("float");
