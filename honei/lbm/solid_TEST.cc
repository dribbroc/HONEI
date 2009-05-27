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

            std::cout << target << std::endl;

            //FluidToSolid test:
            DenseMatrix<bool> target_fts(20, 20, false);
            DenseMatrix<bool> tm1_fts(20, 20, false);
            DenseMatrix<bool> t_fts(20, 20, false);

            Line<DataType_, lbm_solid_dims::D2> line_3(DataType_(1), DataType_(1), DataType_(10), DataType_(10));
            Line<DataType_, lbm_solid_dims::D2> line_4(DataType_(10), DataType_(10), DataType_(15), DataType_(1));
            Line<DataType_, lbm_solid_dims::D2> line_5(DataType_(15), DataType_(1), DataType_(1), DataType_(1));

            Polygon<DataType_, lbm_solid_dims::D2> tri_1(3);
            tri_1.add_line(line_3);
            tri_1.add_line(line_4);
            tri_1.add_line(line_5);

            tri_1.value();

            Line<DataType_, lbm_solid_dims::D2> line_6(DataType_(2), DataType_(2), DataType_(11), DataType_(11));
            Line<DataType_, lbm_solid_dims::D2> line_7(DataType_(11), DataType_(11), DataType_(16), DataType_(2));
            Line<DataType_, lbm_solid_dims::D2> line_8(DataType_(16), DataType_(2), DataType_(2), DataType_(2));

            Polygon<DataType_, lbm_solid_dims::D2> tri_2(3);
            tri_2.add_line(line_6);
            tri_2.add_line(line_7);
            tri_2.add_line(line_8);

            tri_2.value();

            ScanConversion<Tag_>::value(tri_1, tm1_fts, DataType_(1), DataType_(1));
            std::cout << "Matrix at time t-1: " << std::endl;
            std::cout << tm1_fts << std::endl;
            ScanConversion<Tag_>::value(tri_2, t_fts, DataType_(1), DataType_(1));
            std::cout << "Matrix at time t: " << std::endl;
            std::cout << t_fts << std::endl;

            FluidToSolidCells<Tag_>::value(tm1_fts, t_fts, target_fts);
            std::cout << "Fluid to solid cells: " << std::endl;
            std::cout << target_fts << std::endl;

            //Solid to fluid test:
            DenseMatrix<bool> target_stf(20, 20, false);
            SolidToFluidCells<Tag_>::value(tm1_fts, t_fts, target_stf);
            std::cout << "Solid to fluid cells: " << std::endl;
            std::cout << target_stf << std::endl;

            //STFExtrapolation test:
            DenseMatrix<DataType_> h(20, 20, DataType_(5));
            DenseMatrix<DataType_> u(20, 20, DataType_(-5));
            DenseMatrix<DataType_> v(20, 20, DataType_(-5));

            STFExtrapolation<Tag_, lbm_solid_extrapolation_methods::SIMPLE>::value(t_fts, target_stf, h, u, v, DataType_(1), DataType_(1), DataType_(1), DataType_(1));

            std::cout << "h after STF extrapolation: " << std::endl;
            std::cout << h << std::endl;
            std::cout << "u after STF extrapolation: " << std::endl;
            std::cout << u << std::endl;
            std::cout << "v after STF extrapolation: " << std::endl;
            std::cout << v << std::endl;

            //FTSExtrapolation test:
            DenseMatrix<DataType_> h_2(20, 20, DataType_(5));
            DenseMatrix<DataType_> u_2(20, 20, DataType_(-5));
            DenseMatrix<DataType_> v_2(20, 20, DataType_(-5));

            FTSExtrapolation<Tag_, lbm_solid_extrapolation_methods::SIMPLE>::value(t_fts, target_fts, h_2, u_2, v_2, DataType_(1), DataType_(1), DataType_(1), DataType_(1));

            std::cout << "h after FTS extrapolation: " << std::endl;
            std::cout << h_2 << std::endl;
            std::cout << "u after FTS extrapolation: " << std::endl;
            std::cout << u_2 << std::endl;
            std::cout << "v after FTS extrapolation: " << std::endl;
            std::cout << v_2 << std::endl;

        }
};
SolidTest<tags::CPU, float> solidtest_float("float");
