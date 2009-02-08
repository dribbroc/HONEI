/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/util/type_traits.hh>
#include <unittest/unittest.hh>

#include <limits>

using namespace honei;
using namespace tests;

template <typename DataType_>
class TypeTraitsTest :
    public QuickTest
{
    public:
        TypeTraitsTest(const std::string & type) :
            QuickTest("type_traits_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned size(99), offset(7), skip(17);

            DataType_ * elements(TypeTraits<DataType_>::allocate(size));

            TEST_CHECK_EQUAL(reinterpret_cast<unsigned long long>(elements) % 16, 0ul);

            TypeTraits<DataType_>::create(elements, size);
            TypeTraits<DataType_>::fill(elements, size, DataType_(7));

            for (DataType_ * i(elements), * i_end(elements + size) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(7), std::numeric_limits<DataType_>::epsilon());
            }

            TypeTraits<DataType_>::fill(elements + offset, size - skip - offset, DataType_(3));

            for (DataType_ * i(elements), * i_end(elements + offset) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(7), std::numeric_limits<DataType_>::epsilon());
            }

            for (DataType_ * i(elements + offset), * i_end(elements + size - skip) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(3), std::numeric_limits<DataType_>::epsilon());
            }

            for (DataType_ * i(elements + size - skip), * i_end(elements + size) ; i != i_end ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, DataType_(7), std::numeric_limits<DataType_>::epsilon());
            }

            TypeTraits<DataType_>::free(elements, size);
            TEST_CHECK(true);

            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(1), 16ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(16), 16ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(17), 32ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(1024), 1024ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(1025), 1040ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(1026), 1040ul);
            TEST_CHECK_EQUAL(intern::multiple_of_sixteen(1025*1025), 1050640ul);
        }
};
TypeTraitsTest<float> type_traits_test_float("float");
TypeTraitsTest<double> type_traits_test_double("double");
