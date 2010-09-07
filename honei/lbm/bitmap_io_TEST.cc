/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
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

#include <honei/lbm/bitmap_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template<typename DT_, typename FileType_>
class BitmapIOTest:
    public BaseTest
{
    private:
        std::string _tag;
    public:
        BitmapIOTest(const std::string & tag) :
            BaseTest("Bitmap I/O test " + tag)
        {
            _tag = tag;
        }

        virtual void run() const
        {
            std::string filename;
            if(_tag.find("PPM") == -1)
                filename = "test_2.pgm";
            else
                filename = "test_2.ppm";

            DT_ scale(1);
            DenseMatrix<DT_> result(BitmapIO<FileType_>::read_scalar_field(filename, scale));

            TEST_CHECK_EQUAL(result[50][50], DT_(0));
            TEST_CHECK_EQUAL(result[0][0], DT_(1));
            TEST_CHECK_EQUAL(result[99][99], DT_(1));
            TEST_CHECK_EQUAL(result[0][99], DT_(1));
            TEST_CHECK_EQUAL(result[99][0], DT_(1));
            //std::cout << result << std::endl;

            std::string outname(_tag);
            outname += ".ppm";
            BitmapIO<FileType_>::write_scalar_field(result, outname);

            if(_tag.find("PPM") != -1)
            {
                DenseMatrix<DT_> result_2(BitmapIO<FileType_>::read_scalar_field(outname, scale));
                TEST_CHECK_EQUAL(result_2[50][50], DT_(0));
                TEST_CHECK_EQUAL(result_2[0][0], DT_(1));
                TEST_CHECK_EQUAL(result_2[99][99], DT_(1));
                TEST_CHECK_EQUAL(result_2[0][99], DT_(1));
                TEST_CHECK_EQUAL(result_2[99][0], DT_(1));
            }
        }
};
BitmapIOTest<float, io_formats::PGM> bitmapio_test_float("float_PGM");
BitmapIOTest<double, io_formats::PGM> bitmapio_test_double("double_PGM");
BitmapIOTest<float, io_formats::PPM> bitmapio_test_float_ppm("float_PPM");
BitmapIOTest<double, io_formats::PPM> bitmapio_test_double_ppm("double_PPM");
