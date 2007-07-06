/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
 
///< \todo Fix compile errors
#include <libla/banded_matrix.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class BandedMatrixCreationTest :
    public BaseTest
{
    public:
        BandedMatrixCreationTest(const std::string & type) :
            BaseTest("banded_matrix_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                std::tr1::shared_ptr<BandedMatrix<DataType_> > bm1(new BandedMatrix<DataType_>(size));
                TEST_CHECK(true);
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(1)));
                std::tr1::shared_ptr<BandedMatrix<DataType_> > bm2(new BandedMatrix<DataType_>(size,dv1));
                TEST_CHECK(true);
                
            }
            
            DenseVector<DataType_> * dv2(new DenseVector<DataType_>(5, static_cast<DataType_>(1)));                
            TEST_CHECK_THROWS (std::tr1::shared_ptr<BandedMatrix<DataType_> > bm3(new BandedMatrix<DataType_>(6, dv2)),
                VectorSizeDoesNotMatch);
        }
};

BandedMatrixCreationTest<float> banded_matrix_creation_test_float("float");
BandedMatrixCreationTest<double> banded_matrix_creation_test_double("double");

template <typename DataType_>
class BandedMatrixQuickTest :
    public QuickTest
{
    public:
        BandedMatrixQuickTest(const std::string & type) :
            QuickTest("banded_matrix_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            std::tr1::shared_ptr<BandedMatrix<DataType_> > bm1(new BandedMatrix<DataType_>(size));
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(1)));
            std::tr1::shared_ptr<BandedMatrix<DataType_> > bm2(new BandedMatrix<DataType_>(size,dv1));
            DenseVector<DataType_> dv3 = bm2->band(1);
            dv3[size/2] = static_cast<DataType_>(5);
            DenseVector<DataType_> * dv4 (new DenseVector<DataType_>(size, static_cast<DataType_>(0)));
            (*dv4)[size/2] = static_cast<DataType_>(5);
            TEST_CHECK_EQUAL(&bm2->band(0), dv1);
            TEST_CHECK_EQUAL(bm2->band(1), *dv4);            
            TEST_CHECK_EQUAL(bm2->rows(), size);
            TEST_CHECK_EQUAL(bm2->columns(), size);            
            DenseVector<DataType_> * dv2(new DenseVector<DataType_>(5, static_cast<DataType_>(1)));                
            TEST_CHECK_THROWS (std::tr1::shared_ptr<BandedMatrix<DataType_> > bm3(new BandedMatrix<DataType_>(6, dv2)),
                VectorSizeDoesNotMatch);            
        }
};
BandedMatrixQuickTest<float>  banded_matrix_quick_test_float("float");
BandedMatrixQuickTest<double> banded_matrix_quick_test_double("double");
