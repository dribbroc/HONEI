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

#include <honei/swe/correction_processing.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>

using namespace tests;
using namespace std;
using namespace boundaries;

template <typename Tag_, typename DT1_>
class CorrectionProcessingTest:
    public BaseTest
{
    public:
        CorrectionProcessingTest(const std::string & tag) :
            BaseTest("CorrectionProcessing test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long d_width = 1;
            unsigned long d_height = 1;

            unsigned long entries(3 * ((d_width * d_height) + 4 * (d_width + d_height + 4)));
            DenseVector<DT1_> u(entries, DT1_(10));
            DenseVector<DT1_> v(entries, DT1_(10));
            DenseVector<DT1_> w(entries, DT1_(10));
            DenseVector<DT1_> pu(entries, DT1_(1));
            DenseVector<DT1_> pv(entries, DT1_(1));
            DenseVector<DT1_> pw(entries, DT1_(1));
            DenseMatrix<DT1_> h(entries, entries, DT1_(0));

            DenseVector<DT1_> ref(entries, DT1_(10));
            for(unsigned long i = 0; i < entries; i++)
            {
                if( i > 35 && i < 39)
                {
                    ref[i] = DT1_(5.5);
                }
            }

            CorrectionProcessing<REFLECT, Tag_>::value(pu, pv, pw, u, v, w, 1, 1, h);
#ifndef HONEI_SSE
            TEST_CHECK_EQUAL(u, ref);
            TEST_CHECK_EQUAL(v, ref);
            TEST_CHECK_EQUAL(w, ref);
#endif

#ifdef HONEI_SSE
            //To normlize testcases for SSE
            for(unsigned long i = 0; i < entries; i++)
            {
                if( i > 35 && i < 39)
                {
                }
                else
                {
                    u[i] = DT1_(10);
                    v[i] = DT1_(10);
                    w[i] = DT1_(10);
                }
            }
#endif


            TEST_CHECK_EQUAL(u, ref);
            TEST_CHECK_EQUAL(v, ref);
            TEST_CHECK_EQUAL(w, ref);
        }
};
CorrectionProcessingTest<tags::CPU, float> correction_test_float("float");
CorrectionProcessingTest<tags::CPU, double> correction_test_double("double");
CorrectionProcessingTest<tags::CPU::MultiCore, float> correction_test_float_mc("float MC");
CorrectionProcessingTest<tags::CPU::MultiCore, double> correction_test_double_mc("double MC");

#ifdef HONEI_SSE
CorrectionProcessingTest<tags::CPU::SSE, float> correction_test_float_sse("float SSE");
CorrectionProcessingTest<tags::CPU::SSE, double> correction_test_double_sse("double SSE");
#endif
