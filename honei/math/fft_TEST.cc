/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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


#include <honei/math/fft.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <math.h>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT_>
class FFTTestSparse:
    public BaseTest
{
    public:
        FFTTestSparse(const std::string & tag) :
            BaseTest("FFT test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT_> x(8, 4711);
            for (unsigned long i(0) ; i < x.size() ; ++i)
                x[i] = sin(DT_(i));
            DenseVector<DT_> y = FFT::forward(x);
            DenseVector<DT_> z = FFT::backward(y);
            TEST_CHECK_EQUAL(z, x);
        }
};
FFTTestSparse<tags::CPU, double> fft_test_double("double");
