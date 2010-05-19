/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de
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

#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/math/ludecomposition.hh>
#include <honei/la/product.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename DT_, typename Tag_>
class LUTest:
    public BaseTest
{
    public:
        LUTest(const std::string & tag) :
            BaseTest("lu_decomp_test_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DT_> a(size, size);
                for (unsigned long i(0) ; i < a.size() ; ++i)
                {
                    a.elements()[i] = DT_(i + 1) / DT_(3);
                }
                DenseVector<DT_> x(size);
                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    x.elements()[i] = DT_(i + 1) / DT_(1);
                }

                DenseVector<DT_> b(Product<Tag_>::value(a, x));

                DenseVector<DT_> result(size);
                std::cout<<b;
                LUDecomposition<Tag_>::value(a, b, result);
                TEST_CHECK_EQUAL(result, x);
            }
        }
};
//LUTest<float, tags::CPU> lu_test_float("float");
//LUTest<double, tags::CPU> lu_test_double("double");

template<typename DT_, typename Tag_>
class LUQuickTest:
    public QuickTest
{
    public:
        LUQuickTest(const std::string & tag) :
            QuickTest("lu_decomp_quick_test_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT_> a(3, 3, DT_(1));
            DenseVector<DT_> b(3, DT_(1));
            a[0][0] = DT_(7);
            a[0][1] = DT_(-2);
            a[0][2] = DT_(0);
            a[1][0] = DT_(-2);
            a[1][1] = DT_(6);
            a[1][2] = DT_(2);
            a[2][0] = DT_(0);
            a[2][1] = DT_(2);
            a[2][2] = DT_(5);

            b[0] = DT_(3);
            b[1] = DT_(3);
            b[2] = DT_(0);


            DenseVector<DT_> x_analytical(3);
            x_analytical[0] = DT_(2./3.);
            x_analytical[1] = DT_(5./6.);
            x_analytical[2] = DT_(-1./3.);

            DenseVector<DT_> result(3);
            LUDecomposition<Tag_>::value(a, b, result);
            TEST_CHECK_EQUAL(result, x_analytical);
        }
};
LUQuickTest<float, tags::CPU> lu_quick_test_float("float");
LUQuickTest<double, tags::CPU> lu_quick_test_double("double");

template<typename DT_, typename Tag_>
class PLUQuickTest:
    public QuickTest
{
    public:
        PLUQuickTest(const std::string & tag) :
            QuickTest("plu_decomp_quick_test_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT_> a(2, 2);
            DenseVector<DT_> b(2);
            a[0][0] = DT_(0);
            a[0][1] = DT_(1);
            a[1][0] = DT_(1);
            a[1][1] = DT_(1);

            b[0] = DT_(1);
            b[1] = DT_(0);


            DenseVector<DT_> x_analytical(2);
            x_analytical[0] = DT_(-1);
            x_analytical[1] = DT_(1);

            DenseVector<DT_> result(2);
            LUDecomposition<Tag_>::value(a, b, result);
            TEST_CHECK_EQUAL(result, x_analytical);
        }
};
//PLUQuickTest<float, tags::CPU> plu_quick_test_float("float");
//PLUQuickTest<double, tags::CPU> plu_quick_test_double("double");
