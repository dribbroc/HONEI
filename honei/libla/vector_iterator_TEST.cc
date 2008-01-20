/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/libla/vector_iterator.hh>
#include <honei/libla/banded_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using namespace tests;

template <typename DataType_>
class BandedMatrixElementIterationTest :
    public BaseTest
{
    public:
        BandedMatrixElementIterationTest(const std::string & type) :
            BaseTest("banded_matrix_vector_iteration_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                BandedMatrix<DataType_> bm(size);

                typename BandedMatrix<DataType_>::VectorIterator b(bm.begin_bands()), b_end(bm.end_bands());
                for (unsigned long i(1 - size) ; i < (size - 1) ; ++i, ++b)
                {
                    TEST_CHECK_EQUAL(b.index(), i);

                    typename BandedMatrix<DataType_>::ConstVectorIterator cb(b);
                    TEST_CHECK_EQUAL(std::distance((*cb).begin_elements(), (*cb).end_elements()), size);

                    for (typename Vector<DataType_>::ConstElementIterator ce((*cb).begin_elements()), ce_end((*cb).end_elements()) ;
                            ce != ce_end ; ++ce)
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*ce, 0, std::numeric_limits<DataType_>::epsilon());
                    }
                }
            }
        }
};

BandedMatrixElementIterationTest<float> banded_matrix_element_iteration_test_float("float");
BandedMatrixElementIterationTest<double> banded_matrix_element_iteration_test_double("double");
