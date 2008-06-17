/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/banded_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using  namespace tests;

template <typename DataType_>
class BandedMatrixQ1QuickTest :
    public QuickTest
{
    public:
        BandedMatrixQ1QuickTest(const std::string & type) :
            QuickTest("banded_matrix_q1_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (1089);
            DenseVector<DataType_> dv0(size, DataType_(4711));
            for (unsigned long index(0) ; index < size ; ++index)
            {
                dv0[index] = DataType_(index);
            }
            DenseVector<DataType_> dv1(size, DataType_(1));
            DenseVector<DataType_> dv2(size, DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(3));
            DenseVector<DataType_> dv4(size, DataType_(4));
            DenseVector<DataType_> dv5(size, DataType_(5));
            DenseVector<DataType_> dv6(size, DataType_(6));
            DenseVector<DataType_> dv7(size, DataType_(7));
            DenseVector<DataType_> dv8(size, DataType_(8));
            DenseVector<DataType_> dv9(size, DataType_(9));

            BandedMatrixQ1<DataType_> bm1(size, dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8, dv9);
            BandedMatrixQ1<DataType_> bm0(size, dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy(), dv0.copy());
            TEST_CHECK_NOT_EQUAL(bm0, bm1);
            TEST_CHECK_EQUAL(bm0, bm0);
            TEST_CHECK_EQUAL(bm1, bm1.copy());

            BandedMatrix<DataType_> am(bm1);
            BandedMatrixQ1<DataType_> bm3(am);
            /// \todo activate when bm==bmq1 is available
            /*TEST_CHECK_EQUAL(am, bm1);
            TEST_CHECK_EQUAL(bm3, am) */;
            TEST_CHECK_EQUAL(bm1, bm3);


            DenseVector<DataType_> dv01(5, DataType_(1));
            TEST_CHECK_THROWS(BandedMatrixQ1<DataType_> bm01(9, dv01, dv01, dv01, dv01, dv01, dv01, dv01, dv01, dv01), VectorSizeDoesNotMatch);
        }
};

BandedMatrixQ1QuickTest<float> banded_matrix_q1_quick_test_float("float");
BandedMatrixQ1QuickTest<double> banded_matrix_q1quick_test_double("double");

