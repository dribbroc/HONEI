/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 - 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/prolongation_matrix.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/la/dense_matrix.hh>
#include <cmath>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class ProlongationMatrixTest:
    public BaseTest
{
    public:
        ProlongationMatrixTest(const std::string & tag) :
            BaseTest("Prolongation matrix assembly test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            SparseMatrix<DT1_> prolmat(25, 9);
            DenseVector<unsigned long> indices_fine(25);
            DenseVector<unsigned long> indices_coarse(9);
            for(unsigned long i(0) ; i < 25 ; ++i)
                indices_fine[i] = i;
            for(unsigned long i(0) ; i < 9 ; ++i)
                indices_coarse[i] = i;

            ProlongationMatrix<Tag_>::value(prolmat, indices_fine, indices_coarse);
            std::cout << prolmat;

            TEST_CHECK_EQUAL(prolmat(0,0), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(1,0), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(1,1), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(2,1), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(3,1), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(3,2), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(4,2), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(5,0), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(5,3), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(6,0), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(6,1), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(6,3), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(6,4), DT1_(0.25));

            TEST_CHECK_EQUAL(prolmat(7,1), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(7,4), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(8,1), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(8,2), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(8,4), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(8,5), DT1_(0.25));

            TEST_CHECK_EQUAL(prolmat(9,2), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(9,5), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(10,3), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(11,3), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(11,4), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(12,4), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(13,4), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(13,5), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(14,5), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(15,3), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(15,6), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(16,3), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(16,4), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(16,6), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(16,7), DT1_(0.25));

            TEST_CHECK_EQUAL(prolmat(17,4), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(17,7), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(18,4), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(18,5), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(18,7), DT1_(0.25));
            TEST_CHECK_EQUAL(prolmat(18,8), DT1_(0.25));

            TEST_CHECK_EQUAL(prolmat(19,5), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(19,8), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(20,6), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(21,6), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(21,7), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(22,7), DT1_(1));

            TEST_CHECK_EQUAL(prolmat(23,7), DT1_(0.5));
            TEST_CHECK_EQUAL(prolmat(23,8), DT1_(0.5));

            TEST_CHECK_EQUAL(prolmat(24,8), DT1_(1));

        }
};
ProlongationMatrixTest<tags::CPU, float> cpu_prolongation_matrix_test_float("float");

