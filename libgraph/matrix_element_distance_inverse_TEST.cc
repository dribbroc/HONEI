/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the Graph C++ library. LibGraph is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libla/dense_matrix.hh>
#include <libgraph/matrix_element_distance_inverse.hh>
#include <libgraph/matrix_element_distance.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using  namespace tests;
template <typename DataType_>
class MatrixElementDistanceInverseQuickTest :
    public QuickTest
{
    public:
        MatrixElementDistanceInverseQuickTest(const std::string & type) :
            QuickTest("matrix_element_distance_inverse_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long columns(3), rows(2);
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm(new DenseMatrix<DataType_>
                (3, 2));

            for (typename MutableMatrix<DataType_>::ElementIterator i(dm->begin_elements()),
                    i_end(dm->end_elements()) ; i != i_end ; ++i)
            {
                *i = i.index();
            }
            DenseMatrix<DataType_> distance = MatrixElementDistance<DataType_>::value(*dm);            
            DenseMatrix<DataType_> distance_inverse = MatrixElementDistanceInverse<DataType_>::value(*dm);

            TEST_CHECK_EQUAL(distance_inverse[0][0], 0);
            TEST_CHECK_EQUAL(distance_inverse[1][1], 0);
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[0][1], 1 / distance[0][1], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[1][0], 1 / distance[1][0], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[0][2], 1 / distance[0][2], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[1][2], 1 / distance[1][2], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[2][0], 1 / distance[2][0], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[2][1], 1 / distance[2][1], 
                std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance_inverse[2][2], 0, 
                std::numeric_limits<DataType_>::epsilon());

            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm2(new DenseMatrix<DataType_>
                (2, 3));
            TEST_CHECK_THROWS(MatrixElementDistanceInverse<DataType_>::value(*dm2), MatrixRowsDoNotMatch);

        }
};
MatrixElementDistanceInverseQuickTest<float>  matrix_element_distance_inverse_quick_test_float("float");
MatrixElementDistanceInverseQuickTest<double> matrix_element_distance_inverse_quick_test_double("double");
