/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#include <honei/libswe/post_processing.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <iostream>

using namespace tests;
using namespace std;
using namespace output_types;

template <typename Tag_, typename DT1_>
class PostProcessingTest:
    public BaseTest
{
    public:
        PostProcessingTest(const std::string & tag) :
            BaseTest("PostProcessing test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT1_> height(100, 100, DT1_(1.2345));
            PostProcessing<output_types::GNUPLOT>::value(height, 1, 100, 100, 1);
            remove("out1.dat");
            TEST_CHECK(true);
        }
};
PostProcessingTest<tags::CPU, float> source_test_float("float");
PostProcessingTest<tags::CPU, double> source_test_double("double");

