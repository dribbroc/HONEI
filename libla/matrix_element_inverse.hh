/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_ELEMENT_INVERSE_HH
#define LIBLA_GUARD_MATRIX_ELEMENT_INVERSE_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>

/**
 * \file
 *
 * Templatized definitions of inverting a matrix's elements.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{
    /**
     * VectorAbsolute is the class template that inverts a matrix's elements.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixElementInverse
    {
        /**
         * Return a matrix's inverse elements. All elements that equal zero will
         * be invariant under this operation.
         *
         * \param matrix DenseMatrix whose non-zero elements shall be inverted.
         **/
        template <typename DataType_> static DenseMatrix<DataType_> & value(DenseMatrix<DataType_> & matrix)
        {
            for (typename MutableMatrix<DataType_>::ElementIterator i(matrix.begin_elements()), i_end(matrix.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i == static_cast<DataType_>(0))
                    continue;

                *i = DataType_(1) / *i;
            }

            return matrix;
        }

		/**
         * Return a matrix's inverse elements. All elements that equal zero will
         * be invariant under this operation.
         *
         * \param matrix SparseMatrix whose non-zero elements shall be inverted.
         **/
        template <typename DataType_> static SparseMatrix<DataType_> & value(SparseMatrix<DataType_> & matrix)
        {
            for (typename MutableMatrix<DataType_>::ElementIterator i(matrix.begin_non_zero_elements()), i_end(matrix.end_non_zero_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(1) / *i;
            }

            return matrix;
        }

        /**
         * Return a matrix's inverse elements. All elements that equal zero will
         * be invariant under this operation.
         *
         * \param matrix BandedMatrix whose non-zero elements shall be inverted.
         **/
		template <typename DataType_> static BandedMatrix<DataType_> & value(BandedMatrix<DataType_> & matrix)
        {
			for (typename BandedMatrix<DataType_>::VectorIterator i(matrix.begin_bands()), 
					i_end(matrix.end_bands()) ; i != i_end ; ++i)
			{
				DenseVector<DataType_> band = *i;
				int middle_index = matrix.rows() -1;
				// If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                if (i.index() >= middle_index)
				{
					unsigned long end = band.size() - (i.index() - middle_index); //Calculation of the element-index to stop in iteration!
					for(typename Vector<DataType_>::ElementIterator b(band.begin_elements()), b_end(band.element_at(end)) ; b != b_end ; ++b)
					{
						if (*b == static_cast<DataType_>(0))
							continue;

						*b = DataType_(1) / *b;
					}
				}
				else
				{
					unsigned long start = middle_index - i.index(); //Calculation of the element-index to start in iteration!
					for(typename Vector<DataType_>::ElementIterator b(band.element_at(start)), b_end(band.end_elements()) ; b != b_end ; ++b)
                    {
						if (*b == static_cast<DataType_>(0))
							continue;

						*b = DataType_(1) / *b;
					}
				}
			}

            return matrix;
        }

    };
}

#endif
