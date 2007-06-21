/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach  <sven.mallach@uni-dortmund.de>
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
#include <libgraph/matrix_element_distance.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;
template <typename DataType_>
class MatrixElementDistanceQuickTest :
    public QuickTest
{
    public:
        MatrixElementDistanceQuickTest(const std::string & type) :
            QuickTest("matrix_element_distance_quick_test<" + type + ">")
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
            DenseMatrix<DataType_> distance = MatrixElementDistance::value(*dm);
            
            TEST_CHECK_EQUAL(distance[0][0], 0);
            TEST_CHECK_EQUAL(distance[1][1], 0);            
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][1], (8), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][0], (8), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[0][2], (32), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(distance[1][2], (8), std::numeric_limits<DataType_>::epsilon());  
            
            std::tr1::shared_ptr<DenseMatrix<DataType_> > dm2(new DenseMatrix<DataType_>
                (3, 3));            
            TEST_CHECK_THROWS(MatrixElementDistance::value(*dm2), MatrixRowsDoNotMatch);
                                                                                              
        }
};
MatrixElementDistanceQuickTest<float>  matrix_element_distance_quick_test_float("float");
MatrixElementDistanceQuickTest<double> matrix_element_distance_test_quick_double("double");
