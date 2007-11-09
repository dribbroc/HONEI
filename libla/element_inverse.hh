/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_ELEMENT_INVERSE_HH
#define LIBLA_GUARD_ELEMENT_INVERSE_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

namespace honei
{
    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_ = tags::CPU> struct ElementInverse
    {
        /**
         * \name Element inversions
         * \{
         *
         * \brief Returns the inverse values of all of an entity's elements.
         *
         * \param x The entity whose elements' inverse values shall be computed.
         *
         * \retval x Will modify the entity x and return it.
         */

        template <typename DataType_>
        static DenseVector<DataType_> & value(DenseVector<DataType_> & x)
        {
            CONTEXT("When calculating the inverse DenseVector elements");

            for (typename Vector<DataType_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                 if (*l == static_cast<DataType_>(0))
                    continue;

                *l = DataType_(1) / *l;
            }
            return x;
        }

        template <typename DataType_>
        static SparseVector<DataType_> & value(SparseVector<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseVector elements");


            for (typename Vector<DataType_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l = DataType_(1) / *l;
            }
            return x;
        }

        template <typename DataType_>
        static DenseMatrix<DataType_> & value(DenseMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseMatrix elements");

            for (typename MutableMatrix<DataType_>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i == static_cast<DataType_>(0))
                    continue;

                *i = DataType_(1) / *i;
            }

            return x;
        }

        template <typename DataType_>
        static SparseMatrix<DataType_> & value(SparseMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseMatrix elements");

            for (typename MutableMatrix<DataType_>::ElementIterator i(x.begin_non_zero_elements()),
                    i_end(x.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                *i = DataType_(1) / *i;
            }

            return x;
        }

        template <typename DataType_>
        static BandedMatrix<DataType_> & value(BandedMatrix<DataType_> & x)
        {
            CONTEXT("When calculating the inverse value of BandedVector elements");

            for (typename BandedMatrix<DataType_>::VectorIterator i(x.begin_bands()),
                    i_end(x.end_bands()) ; i != i_end ; ++i)
            {
                DenseVector<DataType_> band = *i;
                int middle_index = x.rows() -1;
                // If we are above or on the diagonal band, we start at Element 0 and go on
                // until Element band_size-band_index.
                if (i.index() >= middle_index)
                {
                    //Calculation of the element-index to stop in iteration!
                    unsigned long end = band.size() - (i.index() - middle_index);

                    for (typename Vector<DataType_>::ElementIterator b(band.begin_elements()),
                            b_end(band.element_at(end)) ; b != b_end ; ++b)
                    {
                        if (*b == DataType_(0))
                            continue;

                        *b = DataType_(1) / *b;
                    }
                }
                else
                {
                    //Calculation of the element-index to start in iteration!
                    unsigned long start = middle_index - i.index();
                    for (typename Vector<DataType_>::ElementIterator b(band.element_at(start)),
                            b_end(band.end_elements()) ; b != b_end ; ++b)
                    {
                        if (*b == DataType_(0))
                            continue;

                        *b = DataType_(1) / *b;
                    }
                }
            }

            return x;
        }

        /// \}

        #ifdef BENCHM
        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(unsigned long a_rows, unsigned long a_columns = 1, double nonzero = 1)
        {
            BenchmarkInfo result;
            result.flops = 0;
            result.load = 0;
            result.store = 0;
            cout << endl << "!! No detailed benchmark info available !!" << endl;

            return result; 
        }
        #endif
    };
}
#endif
