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

#include <libswe/flow_processing.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

using namespace tests;
using namespace std;
using namespace directions;
template <typename Tag_, typename DT1_>
class FlowProcessingTest:
    public BaseTest
{
    public:
        FlowProcessingTest(const std::string & tag) :
            BaseTest("FlowProcessing test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> vector_x(12000, DT1_(1));
            DenseVector<DT1_> vector_y(12000, DT1_(1));
            DenseVector<DT1_> result_x = FlowProcessing<tags::CPU, X>::value(vector_x);
            FlowProcessing<tags::CPU, Y>::value(vector_y);
            DenseVector<DT1_> analytical_result_x(12000, DT1_(1));
            DenseVector<DT1_> analytical_result_y(12000, DT1_(1));
            for(unsigned long i = 0; i < 12000; ++i)
            {
                if((i + 2) % 3 == 0)
                    analytical_result_x[i] = DT1_(5.905);

                if((i + 1) % 3 == 0)
                    analytical_result_y[i] = DT1_(5.905);
            }

            TEST_CHECK_EQUAL(vector_x, analytical_result_x);
            TEST_CHECK_EQUAL(vector_y, analytical_result_y);
        }
};
FlowProcessingTest<tags::CPU, float> flow_test_float("float");
FlowProcessingTest<tags::CPU, double> flow_test_double("double");

