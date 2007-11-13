/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libmath/parser.hh>
#include <unittest/unittest.hh>

using namespace tests;
using namespace std;

template <typename DT_>
class ParserTest:
    public BaseTest
{
    public:
        ParserTest(const std::string & tag) :
            BaseTest("Parser Test <" + tag + ">")
        {
            //register_tag(Tag_::name);
        }

        virtual void run() const
        {
            char * t = new char[10];
            //t = "123.123";
            DT_ f = Parser<DT_>::parse(t);
            //cout << f << endl;

            DenseVector<DT_> result = Parser<DT_>::parse("parser_test_float.txt", 9);
            for(int i = 0; i < 9; i++)
            {
                TEST_CHECK_EQUAL(result[i],i);
            }

            DenseMatrix<DT_> resultm = Parser<DT_>::parse("parser_test_float_m.txt", 9, 3);
            for(int i = 0; i < 9; i++)
            {
                TEST_CHECK_EQUAL(resultm[0][i],i);
            }


            TEST_CHECK(true);
        }
};
ParserTest<float> parser_test_float("float");
ParserTest<double> parser_test_double("double");

