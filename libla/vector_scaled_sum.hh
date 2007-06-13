 /* vim: set sw=4 sts=4 et nofoldenable : */
 
 /*
  * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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
  */

#ifndef LIBLA_GUARD_VECTOR_SCALED_SUM_HH
#define LIBLA_GUARD_VECTOR_SCALED_SUM_HH 1

#include <libla/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Templatized definitions of linearcombinations (in libla notation: VectorScaledSum)<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{
	/**
      * VectorScaledSum is the class template for the scaled vector +/- scaled vector operation.
      *
      * \ingroup grpvectoroperations
      **/
	template <typename Tag_ = tags::CPU>
    struct VectorScaledSum
	{
        /**
         * Returns the the resulting vector of the linear lombination of two given DenseVector instances and two 
         * scalars.
         * \param leftScal The scalar provided in order to be multiplied with left
         * \param rightScal analogous
         **/
        template < typename DtLeft_ , typename DtRight_, typename DtLeftScal_, typename DtRightScal_> 
        static DenseVector<DtLeft_> value(DenseVector<DtLeft_> & left, 
                        const DenseVector<DtRight_> & right,          
                        const DtLeftScal_ & leftscal, const DtRightScal_ & rightscal)
		{
			if (left.size() != right.size())
            	                throw VectorSizeDoesNotMatch(right.size(), left.size());
            
                        //DenseVector<DataType_> result(left.size(), 0, 0, 1); old style
            
                        for (typename Vector<DtLeft_>::ElementIterator l(left.begin_elements()), 
                                l_end(left.end_elements()) ; l != l_end ; ++l)
                        {
                                left[l.index()] = ((*l) * leftscal) + (rightscal * right[l.index()]);
                        }
                        return left;
		}
		        /**
          * Returns the the resulting vector of the lincomb. of two given SparseVector instances and two 
          * scalars.
          * \param leftScal The scalar provided in order to be multiplied with left
          * \param rightScal analogous
          **/
        template < typename DtLeft_ , typename DtRight_, typename DtLeftScal_, typename DtRightScal_>
		static SparseVector<DtLeft_> value(SparseVector<DtLeft_> & left, 
                        const SparseVector<DtRight_> & right, const DtLeftScal_ & leftscal, const DtRightScal_ & rightscal)
		{
			if (left.size() != right.size())
                    	        throw VectorSizeDoesNotMatch(right.size(), left.size());
            
                        //SparseVector<DataType_> result(left.size(), right.used_elements() + left.used_elements()); old style
            
                        for (typename Vector<DtLeft_>::ConstElementIterator l(left.begin_non_zero_elements()),
                                l_end(left.end_non_zero_elements()), r(right.begin_non_zero_elements()),
                                r_end(right.end_non_zero_elements()); l != l_end ; ++l)
                        {
            	                if (r.index() < l.index())
                                {
                                        left[r.index()] = (*r) * rightscal;
                                        ++r;
                                }
                                else if (l.index() < r.index())
                                {
                                        left[l.index()] = (*l) * leftscal;
                                        ++l;
                                }
                                else
                                {
                                        left[l.index()] = ((*l) * leftscal) + ((*r) * rightscal);
                                        ++l; 
                                        ++r;
                                }
                                return left;
                         }
		}
		
        /**
          * Returns the the resulting vector of the lincomb. of two given Vector instances (one dense, the other
          * sparse)and two 
          * scalars.
          * \param leftscal The scalar provided in order to be multiplied with left
          * \param rightscal see leftscal
          **/
        template < typename DtLeft_ , typename DtRight_, typename DtLeftScal_, typename DtRightScal_>
		static DenseVector<DtLeft_> value(const DenseVector<DtLeft_> & left, SparseVector<DtRight_> right,
                        const DtRightScal_ & leftscal, const DtRightScal_ & rightscal)
                {
                        if (left.size() != right.size())
                                throw VectorSizeDoesNotMatch(right.size(), left.size());

                        //DenseVector<DataType_> result(left.size(),0, 0, 1); old style

                        for (typename Vector<DataType_>::ConstElementIterator l(left.begin_elements()),
                                        l_end(left.end_elements()), r(right.begin_non_zero_elements()),
                                        r_end(right.end_non_zero_elements()) ;
                                        r != r_end ; )
                        {
                                while (l.index() < r.index() && (l != l_end))
                                {
                                        left[l.index()] = (*l) * leftscal;
                                        ++l;
                                }

                                left[l.index()] = ((*l) * leftscal) + ((*r) * rightscal);
                                ++r;
                        }
                        return left;
           }
	};
	
}

#endif 
