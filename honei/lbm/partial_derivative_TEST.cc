/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/lbm/tags.hh>
#include <unittest/unittest.hh>
#include <honei/lbm/partial_derivative.hh>
#include <iostream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>

class PartialDerivativeTest :
    public TaggedTest<Tag_>
{
    public:
        PartialDerivativeTest(const std::string & type) :
            TaggedTest<Tag_>("partial_derivative_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long height(10), width(10);
            DenseMatrix<DataType_> f(height, width, DataType_(1));
            f(5,5) = 10.;
            DataType_ h(1);

            DenseMatrix<DataType_> result(PartialDerivative<Tag_, X, CENTRALDIFF>::value(f , h));

            std::cout << result << std::endl;

            DenseMatrix<DataType_> result_2(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(f , h));

            std::cout << result_2 << std::endl;
            TEST_CHECK(true);
        }
};
PartialDerivativeTest<tags::CPU, double> pd_test_double("double");
