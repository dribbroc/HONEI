/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/trace.hh>
#include <honei/la/dense_matrix.hh>
#include <unittest/unittest.hh>

#include <limits>

using namespace honei;
using namespace tests;

template <typename DT_>
class DenseMatrixTraceValueTest :
    public BaseTest
{
    public:
        DenseMatrixTraceValueTest(const std::string & type) :
            BaseTest("dense_matrix_trace_value_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DT_ result(0);

                unsigned long rows(size), columns(size + 3);
                DenseMatrix<DT_> dm(rows, columns);

                for (unsigned long i(0), i_end(rows) ; i != i_end ; ++i)
                {
                    for (unsigned long j(0), j_end(columns) ; j != j_end ; ++j)
                    {
                        dm(i, j) = DT_(1.23456789) * (i + 1) * (j + 1);

                        if (i == j)
                            result += DT_(1.23456789) * (i + 1) * (j + 1);
                    }
                }

                DT_ trace(Trace<>::value(dm));
                TEST_CHECK_EQUAL_WITHIN_EPS(trace, result, std::numeric_limits<DT_>::epsilon() * result);
            }
        }
};

DenseMatrixTraceValueTest<float> dense_matrix_trace_value_test_float("float");
DenseMatrixTraceValueTest<double> dense_matrix_trace_value_test_double("double");

