/* vim: set sw=4 sts=4 et nofoldenable : */

/*
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

#ifndef LIBLA_GUARD_ELEMENT_PRODUCT_HH
#define LIBLA_GUARD_ELEMENT_PRODUCT_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

namespace honei
{
    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_ = tags::CPU> struct ElementProduct
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectors elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ;)
            {
                if (r.index() < l.index())
                {
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }
                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            return a;
            ///\todo: perhaps sparsify - in case l.index < r.index Write of 0 possible.
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVector elements");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                *l *= b[l.index()];
            }
            return a;
            ///\todo: perhaps sparsify - if *b[l.index()] == 0 -> write of zero.
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and DenseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()); l != l_end ; ++l)
            {
                *l *= b[l.row()][l.column()];
            }

            return a; ///\todo: perhaps sparsify, dense_matrix[row][col] may be zero.
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                else if (r.index() < l.index())
                {
                    ++r;
                }

                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - in case l.index < r.index set to zero possible.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector elements");

            if (a.rows() != b.rows())
            {
                throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands());
            for (typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()),
                    l_end(a.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists()) //won't check l.exists() here, cause it is always created by Iterator.
                {
                    ++r;
                    continue;
                }

                *l = ElementProduct<>::value(*l, *r);
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a DenseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a SparseMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()); r != r_end ; )
            {
                while (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                *l *= *r;
                ++l; ++r;
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l_end(a.end_elements()) ;
                l != l_end ; ++l)
            {
                *l = DT1_(0);
            }
            /// \todo Left is complete dense at this point.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix and BandedMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and BandedMatrix elements");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                while (r.index() < l.index())
                    ++r;

                *l *= *r;
                ++r; ++l;
            }

            return a;
        }

        /// \}

        #ifdef BENCHM
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(unsigned long rows, unsigned long columns = 1, double nonzero_a = 1, double nonzero_b = 1)
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

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::Cell>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);
        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        /// \}
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::CPU::SSE>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & a, const DenseVector<double> & b);

        /// \}
    };
}
#endif
