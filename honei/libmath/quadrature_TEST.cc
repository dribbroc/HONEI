/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/libmath/quadrature.hh>
#include <unittest/unittest.hh>

#include <cmath>
#include <limits>
#include <iostream>

using namespace honei;
using namespace tests;

namespace
{
    struct Exp :
        public std::unary_function<float, float>
    {
        Exp()
        {
        }

        inline float operator() (const float & x) const
        {
            return ::exp(x);
        }
    };
}

template <typename DataType_, typename QTag_>
class QuadratureTest :
    public BaseTest
{
    public:
        QuadratureTest(const std::string & type, const std::string & tag) :
            BaseTest("quadrature_test<" + type + ", " + tag + ">")
        {
        }

        virtual void run() const
        {
            DataType_ step_size(0.05);
            // Integral of exp over [0,1] with step size 0.1 should be e - 1
            DataType_ v(make_quadrature<tags::Trapezoid>(Exp(), step_size)(DataType_(0), DataType_(1))),
                  s(::exp(1) - ::exp(0));

            TEST_CHECK_EQUAL_WITHIN_EPS(v, s, 1 / DataType_(12 * 12 * 12));
        }
};

QuadratureTest<float, tags::Trapezoid> quadrature_test_float("float", "tags::Trapezoid");
QuadratureTest<double, tags::Trapezoid> quadrature_test_double("double", "tags::Trapezoid");
