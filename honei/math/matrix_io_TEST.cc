/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/matrix_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template<typename DT_>
class MMIOTest:
    public BaseTest
{
    public:
        MMIOTest(const std::string & tag) :
            BaseTest("Matrix Market Format I/O test")
        {
        }

        virtual void run() const
        {
            unsigned long c, r, nz, d;
            std::string filename("5pt_10x10.mtx");
            MatrixIO::get_sizes(filename, c, r, nz, d);
            std::cout << c << " " << r << " " << nz << std::endl;

            DenseMatrix<DT_> matrix = MatrixIO::read_matrix(filename, DT_(0));

            std::cout << matrix;

            TEST_CHECK(true);
        }
};
MMIOTest<float> mmio_test_float_dense("float");

