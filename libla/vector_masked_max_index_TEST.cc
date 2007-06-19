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
 
#include <libla/dense_vector.hh>
#include <libla/vector_masked_max_index.hh>
#include <libla/vector_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorMaskedMaxIndexTest :
    public BaseTest
{
    public:
        DenseVectorMaskedMaxIndexTest(const std::string & type) :
            BaseTest("dense_vector_masked_max_index_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                    size, false));
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
                
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), 
                    i_end(dv->end_elements()) ; i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                
                (*dv)[2] = static_cast<DataType_>(size * 10);
                                    
                for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                    i_end(mask->end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 2 == 0) 
                    {
                        *i=true;
                    }
                }
                unsigned long result(VectorMaskedMaxIndex<DataType_>::value(*dv, *mask));
                TEST_CHECK_EQUAL(result, 2);
            }

            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMaskedMaxIndex<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMaxIndexTest<float> dense_vector_masked_max_index_test_float("float");
DenseVectorMaskedMaxIndexTest<double> dense_vector_masked_max_index_test_double("double"); 

template <typename DataType_>
class DenseVectorMaskedMaxIndexQuickTest :
    public QuickTest
{
    public:
        DenseVectorMaskedMaxIndexQuickTest(const std::string & type) :
            QuickTest("dense_vector_masked_max_index_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<bool> > mask(new DenseVector<bool>(
                size, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
            
            for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), 
                i_end(dv->end_elements()) ; i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }
            
            (*dv)[2] = static_cast<DataType_>(size * 10);
                                
            for (typename Vector<bool>::ElementIterator i(mask->begin_elements()), 
                i_end(mask->end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 2 == 0) 
                {
                    *i=true;
                }
            }
            unsigned long result(VectorMaskedMaxIndex<DataType_>::value(*dv, *mask));
            TEST_CHECK_EQUAL(result, 2);


            std::tr1::shared_ptr<Vector<bool> > mask1(new DenseVector<bool>
                (2, false));
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>
                (3, static_cast<DataType_>(1)));

            TEST_CHECK_THROWS(VectorMaskedMaxIndex<DataType_>::value(*dv1, *mask1), VectorSizeDoesNotMatch);
        }
};
DenseVectorMaskedMaxIndexQuickTest<float>  dense_vector_masked_max_index_quick_test_float("float");
DenseVectorMaskedMaxIndexQuickTest<double> dense_vector_masked_max_index_quick_test_double("double");
