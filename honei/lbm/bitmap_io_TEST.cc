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
    public:
        BitmapIOTest(const std::string & tag) :
            BaseTest("Bitmap I/O test")
        {
        }

        virtual void run() const
        {
            std::string filename("test.bmp");
            DT_ scale(1);
            DenseMatrix<DT_> result(BitmapIO<FileType_>::read_scalar_field(filename, scale));
        }
};
BitmapIOTest<float, io_formats::PNM> bitmapio_test_float("float, BMP");
BitmapIOTest<double, io_formats::PNM> bitmapio_test_double("float, BMP");
