/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/math/sqrt.hh>
#include <honei/util/unittest.hh>

#include <cmath>
#include <limits>

using namespace honei;
using namespace tests;

template <typename DT_, typename Tag_>
class SqrtBigValuesTest :
    public QuickTest
{
    public:
        SqrtBigValuesTest(const std::string & type, const std::string & tag) :
            QuickTest("sqrt_big_values_test<" + type + ", " + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            const unsigned size(4096);
            DenseVector<DT_> test(size), reference(size);

            for (typename DenseVector<DT_>::ElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                *i = DT_(1.0) + DT_(1.23456789) * i.index();
                *j = sqrt(*i);
            }

            Sqrt<Tag_>::value(test);

            for (typename DenseVector<DT_>::ConstElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 6 * *j * std::numeric_limits<DT_>::epsilon());
            }
        }
};
SqrtBigValuesTest<float, tags::CPU> sqrt_big_values_test_float("float", "tags::CPU");
SqrtBigValuesTest<double, tags::CPU> sqrt_big_values_test_double("double", "tags::CPU");
#ifdef HONEI_SSE
SqrtBigValuesTest<float, tags::CPU::SSE> sse_sqrt_big_values_test_float("float", "tags::CPU::SSE");
SqrtBigValuesTest<double, tags::CPU::SSE> sse_sqrt_big_values_test_double("double", "tags::CPU::SSE");
#endif
#ifdef HONEI_CELL
SqrtBigValuesTest<float, tags::Cell> cell_sqrt_big_values_test_float("float", "tags::Cell");
#endif

template <typename DT_, typename Tag_>
class SqrtSmallValuesTest :
    public QuickTest
{
    public:
        SqrtSmallValuesTest(const std::string & type, const std::string & tag) :
            QuickTest("sqrt_small_values_test<" + type + ", " + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            const unsigned size(4096);
            DenseVector<DT_> test(size), reference(size);

            for (typename DenseVector<DT_>::ElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                *i = DT_(1) + i.index() * DT_(3) / DT_(size);
                *j = sqrt(*i);
            }

            Sqrt<Tag_>::value(test);

            for (typename DenseVector<DT_>::ConstElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 6 * std::numeric_limits<DT_>::epsilon());
            }
        }
};

SqrtSmallValuesTest<float, tags::CPU> sqrt_small_values_test_float("float", "tags::CPU");
SqrtSmallValuesTest<double, tags::CPU> sqrt_small_values_test_double("double", "tags::CPU");
#ifdef HONEI_SSE
SqrtSmallValuesTest<float, tags::CPU::SSE> sse_sqrt_small_values_test_float("float", "tags::CPU::SSE");
SqrtSmallValuesTest<double, tags::CPU::SSE> sse_sqrt_small_values_test_double("double", "tags::CPU::SSE");
#endif
#ifdef HONEI_CELL
SqrtSmallValuesTest<float, tags::Cell> cell_sqrt_small_values_test_float("float", "tags::Cell");
#endif

template <typename DT_, typename Tag_>
class SqrtTinyValuesTest :
    public QuickTest
{
    public:
        SqrtTinyValuesTest(const std::string & type, const std::string & tag) :
            QuickTest("sqrt_tiny_values_test<" + type + ", " + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            const unsigned size(4096);
            DenseVector<DT_> test(size), reference(size);

            for (typename DenseVector<DT_>::ElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                *i = i.index() * DT_(1) / DT_(size);
                *j = sqrt(*i);
            }

            Sqrt<Tag_>::value(test);

            for (typename DenseVector<DT_>::ConstElementIterator i(test.begin_elements()), i_end(test.end_elements()), j(reference.begin_elements()) ;
                    i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, 2 * std::numeric_limits<DT_>::epsilon());
            }
        }
};

SqrtTinyValuesTest<float, tags::CPU> sqrt_tiny_values_test_float("float", "tags::CPU");
SqrtTinyValuesTest<double, tags::CPU> sqrt_tiny_values_test_double("double", "tags::CPU");
#ifdef HONEI_SSE
SqrtTinyValuesTest<float, tags::CPU::SSE> sse_sqrt_tiny_values_test_float("float", "tags::CPU::SSE");
SqrtTinyValuesTest<double, tags::CPU::SSE> sse_sqrt_tiny_values_test_double("double", "tags::CPU::SSE");
#endif
#ifdef HONEI_CELL
SqrtTinyValuesTest<float, tags::Cell> cell_sqrt_tiny_values_test_float("float", "tags::Cell");
#endif
