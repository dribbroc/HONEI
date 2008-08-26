/* vim: set sw=4 sts=4 et nofoldenable : */

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

#include <honei/la/dense_vector.hh>
#include <honei/la/const_vector.hh>
#include <unittest/unittest.hh>

using namespace honei;
using namespace tests;

template <typename DataType_>
class ConstVectorTest :
    public QuickTest
{
    public:
        ConstVectorTest(const std::string & type) :
            QuickTest("const_vector_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv(size, DataType_(0));
                ConstVector<DataType_> cv(dv);
                TEST_CHECK_EQUAL(cv.size(), dv.size());
                TEST_CHECK_EQUAL(cv.offset(), dv.offset());
                TEST_CHECK_EQUAL(cv[8], dv[8]);
                TEST_CHECK_EQUAL(cv[size - 2], dv[size - 2]);

                ConstVector<DataType_> cv2(cv);
                TEST_CHECK_EQUAL(cv2.size(), cv.size());
                TEST_CHECK_EQUAL(cv2.offset(), cv.offset());
                TEST_CHECK_EQUAL(cv2[8], cv[8]);
                TEST_CHECK_EQUAL(cv2[size - 2], cv[size - 2]);
                TEST_CHECK_EQUAL(cv, cv2);

                typename DenseVector<DataType_>::ElementIterator j(dv.begin_elements());
                for (typename ConstVector<DataType_>::ConstElementIterator i(cv.begin_elements()), i_end(cv.end_elements()) ;
                        i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL(*i, *j);
                }
            }
        }
};
ConstVectorTest<float> const_vector_test_float("float");
ConstVectorTest<double> const_vector_test_double("double");
