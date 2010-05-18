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

#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <cmath>
#include <honei/math/reordering.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class ReorderingTest:
    public BaseTest
{
    public:
        ReorderingTest(const std::string & tag) :
            BaseTest("Reordering test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<unsigned long> indices_fine(25);
            DenseVector<unsigned long> indices_coarse(9);
            for(unsigned long i(0) ; i < 25 ; ++i)
                indices_fine[i] = i;
            for(unsigned long i(0) ; i < 9 ; ++i)
                indices_coarse[i] = i;
            DenseVector<unsigned long> indices_fine_ref(indices_fine.copy());
            DenseVector<unsigned long> indices_coarse_ref(indices_coarse.copy());

            for(unsigned long i(0) ; i < 25 ; ++i)
                indices_fine[i] = 4711;
            for(unsigned long i(0) ; i < 9 ; ++i)
                indices_coarse[i] = 4711;

            Reordering<Tag_, methods::NATURAL>::value(indices_fine, indices_coarse);

            TEST_CHECK_EQUAL(indices_fine, indices_fine_ref);
            TEST_CHECK_EQUAL(indices_coarse, indices_coarse_ref);

            indices_fine_ref[0] = 0;
            indices_fine_ref[1] = 9;
            indices_fine_ref[2] = 1;
            indices_fine_ref[3] = 10;
            indices_fine_ref[4] = 2;
            indices_fine_ref[5] = 11;
            indices_fine_ref[6] = 12;
            indices_fine_ref[7] = 13;
            indices_fine_ref[8] = 14;
            indices_fine_ref[9] = 15;
            indices_fine_ref[10] = 3;
            indices_fine_ref[11] = 16;
            indices_fine_ref[12] = 4;
            indices_fine_ref[13] = 17;
            indices_fine_ref[14] = 5;
            indices_fine_ref[15] = 18;
            indices_fine_ref[16] = 19;
            indices_fine_ref[17] = 20;
            indices_fine_ref[18] = 21;
            indices_fine_ref[19] = 22;
            indices_fine_ref[20] = 6;
            indices_fine_ref[21] = 23;
            indices_fine_ref[22] = 7;
            indices_fine_ref[23] = 24;
            indices_fine_ref[24] = 8;
            Reordering<Tag_, methods::TWO_LEVEL>::value(indices_fine, indices_coarse);
            TEST_CHECK_EQUAL(indices_fine, indices_fine_ref);
            TEST_CHECK_EQUAL(indices_coarse, indices_coarse_ref);
        }
};
ReorderingTest<tags::CPU, float> cpu_prolongation_matrix_test_float("float");

