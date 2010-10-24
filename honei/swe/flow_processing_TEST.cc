/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
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

#include <honei/swe/flow_processing.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
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
            unsigned long size(12003);
            DenseVector<DT1_> vector_x(size, DT1_(1));
            DenseVector<DT1_> vector_y(size, DT1_(1));
            FlowProcessing<X, Tag_>::value(vector_x);
            FlowProcessing<Y, Tag_>::value(vector_y);
            DenseVector<DT1_> analytical_result_x(size, DT1_(1));
            DenseVector<DT1_> analytical_result_y(size, DT1_(1));
            for(unsigned long i = 0; i < size; ++i)
            {
                if((i + 2) % 3 == 0)
                    analytical_result_x[i] = DT1_(5.905);

                if((i + 1) % 3 == 0)
                    analytical_result_y[i] = DT1_(5.905);
            }

            for (typename DenseVector<DT1_>::ConstElementIterator i(vector_x.begin_elements()), i_end(vector_x.end_elements()),
                    j(analytical_result_x.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 2 * *j * std::numeric_limits<DT1_>::epsilon());
            }

            for (typename DenseVector<DT1_>::ConstElementIterator i(vector_y.begin_elements()), i_end(vector_y.end_elements()),
                    j(analytical_result_y.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 2 * *j * std::numeric_limits<DT1_>::epsilon());
            }

            DenseVector<DT1_> vector_x_2(size, DT1_(1.23456));
            DenseVector<DT1_> vector_y_2(size, DT1_(5.99999));
            DenseVector<DT1_> cpu_result_x_2(size, DT1_(1.23456));
            DenseVector<DT1_> cpu_result_y_2(size, DT1_(5.99999));

            for(unsigned long i = 0; i < size; ++i)
            {
                if((i + 2) % 3 == 0)
                {
                    vector_x_2[i] = DT1_(1.987654321);
                    cpu_result_x_2[i] = DT1_(1.987654321);
                }
                if((i + 1) % 3 == 0)
                {
                    vector_y_2[i] = DT1_(0.123456789);
                    cpu_result_y_2[i] = DT1_(0.123456789);
                }
            }
            FlowProcessing<X, Tag_>::value(vector_x_2);
            FlowProcessing<Y, Tag_>::value(vector_y_2);
            FlowProcessing<X, tags::CPU>::value(cpu_result_x_2);
            FlowProcessing<Y, tags::CPU>::value(cpu_result_y_2);

            for (typename DenseVector<DT1_>::ConstElementIterator i(vector_x_2.begin_elements()), i_end(vector_x_2.end_elements()),
                    j(cpu_result_x_2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 2 * *j * std::numeric_limits<DT1_>::epsilon());
            }

            for (typename DenseVector<DT1_>::ConstElementIterator i(vector_y_2.begin_elements()), i_end(vector_y_2.end_elements()),
                    j(cpu_result_y_2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 2 * *j * std::numeric_limits<DT1_>::epsilon());
            }

            DenseVector<DT1_> vector_x_3(size, DT1_(1.23456));
            DenseVector<DT1_> vector_y_3(size, DT1_(5.99999));
            DenseVector<DT1_> cpu_result_x_3(size, DT1_(1.23456));
            DenseVector<DT1_> cpu_result_y_3(size, DT1_(5.99999));

            for(unsigned long i = 0; i < size; ++i)
            {
                if((i + 2) % 3 == 0)
                {
                    vector_x_3[i] = DT1_(- 1.987654321);
                    cpu_result_x_3[i] = DT1_(- 1.987654321);
                }
                if((i + 1) % 3 == 0)
                {
                    vector_y_3[i] = DT1_(1000);
                    cpu_result_y_3[i] = DT1_(1000);
                }
            }
            FlowProcessing<X, Tag_>::value(vector_x_3);
            FlowProcessing<Y, Tag_>::value(vector_y_3);
            FlowProcessing<X, tags::CPU>::value(cpu_result_x_3);
            FlowProcessing<Y, tags::CPU>::value(cpu_result_y_3);
            for(unsigned long i(0); i < size; ++i )
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(vector_x_3[i], cpu_result_x_3[i], 2 * fabs(cpu_result_x_3[i])*std::numeric_limits<DT1_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(vector_y_3[i], cpu_result_y_3[i], 2 * fabs(cpu_result_y_3[i])* std::numeric_limits<DT1_>::epsilon());
            }
        }
};
FlowProcessingTest<tags::CPU, float> flow_test_float("float");
FlowProcessingTest<tags::CPU, double> flow_test_double("double");
FlowProcessingTest<tags::CPU::MultiCore, double> flow_test_double_mc("double MC");
FlowProcessingTest<tags::CPU::MultiCore, float> flow_test_float_mc("float MC");
#ifdef HONEI_CELL
FlowProcessingTest<tags::Cell, float> flow_test_float_cell("float Cell");
#endif
#ifdef HONEI_SSE
FlowProcessingTest<tags::CPU::SSE, float> flow_test_float_sse("float SSE");
FlowProcessingTest<tags::CPU::SSE, double> flow_test_double_sse("double SSE");
FlowProcessingTest<tags::CPU::MultiCore::SSE, float> flow_test_float_mcsse("float MCSSE");
FlowProcessingTest<tags::CPU::MultiCore::SSE, double> flow_test_double_mcsse("double MCSSE");
#endif
