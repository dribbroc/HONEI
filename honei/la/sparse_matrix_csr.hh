/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
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

#pragma once
#ifndef LIBLA_GUARD_SPARSE_MATRIX_CSR_HH
#define LIBLA_GUARD_SPARSE_MATRIX_CSR_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/util/private_implementation_pattern.hh>
#ifdef HONEI_GMP
#include <gmpxx.h>
#endif


namespace honei
{
    // Forward declarations
    template <typename DataType_> class SparseMatrix;
    template <typename DataType_> class SparseMatrixELL;

    /**
     * \brief SparseMatrixCSR is a sparse matrix with its data kept in the CSR format.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class SparseMatrixCSR :
        PrivateImplementationPattern<SparseMatrixCSR<DataType_>, Shared>
    {
        public:
            /// \name Basic operations
            /// \{


            /**
             * Constructor.
             *
             * \param CSR Vectors the matrix will be created from.
             */
            explicit SparseMatrixCSR(const unsigned long rows, const unsigned long columns,
                    const DenseVector<unsigned long> & Aj, const DenseVector<DataType_> & Ax, const DenseVector<unsigned long> & Ar);

            /**
             * Constructor.
             *
             * \param src The SparseMatrix our matrix will be created from.
             */
            explicit SparseMatrixCSR(const SparseMatrix<DataType_> & src);

            /**
             * Constructor.
             *
             * \param src The SparseMatrixELL our matrix will be created from.
             */
            explicit SparseMatrixCSR(const SparseMatrixELL<DataType_> & src);

            /// Copy-constructor.
            SparseMatrixCSR(const SparseMatrixCSR<DataType_> & other);

            /// Destructor.
            ~SparseMatrixCSR();

            /// \}

            /// Returns the number of our columns.
            unsigned long columns() const;

            /// Returns the number of our rows.
            unsigned long rows() const;

            /// Returns our size, equal to rows and columns.
            unsigned long size() const;

            /// Returns out non zero element count.
            unsigned long used_elements() const;

            /// Retrieves our Aj (indices) vector.
            DenseVector<unsigned long> & Aj() const;

            /// Retrieves our Ax (data) vector.
            DenseVector<DataType_> & Ax() const;

            /// Retrieves our Ar (row start) vector.
            DenseVector<unsigned long> & Ar() const;

            /// Retrieves element at (row, column), unassignable.
            const DataType_ operator() (unsigned long row, unsigned long column) const;

            /// Request a memory access lock for our data.
            void lock(LockMode mode) const;

            /// Release a memory access lock for our data.
            void unlock(LockMode mode) const;

            /// Returns a copy of the matrix.
            SparseMatrixCSR copy() const;
    };

    /**
     * Equality operator for SparseMatrixCSR.
     *
     * Compares if corresponding elements of two ell matrices are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const SparseMatrixCSR<DataType_> & a, const SparseMatrixCSR<DataType_> & b);

    /**
     * Output operator for SparseMatrixCSR.
     *
     * Outputs a ell matrix to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSR<DataType_> & matrix);

    extern template class SparseMatrixCSR<float>;

    extern template bool operator== (const SparseMatrixCSR<float> & a, const SparseMatrixCSR<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSR<float> & matrix);

    extern template class SparseMatrixCSR<double>;

    extern template bool operator== (const SparseMatrixCSR<double> & a, const SparseMatrixCSR<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSR<double> & matrix);

#ifdef HONEI_GMP
    extern template class SparseMatrixCSR<mpf_class>;

    extern template bool operator== (const SparseMatrixCSR<mpf_class> & a, const SparseMatrixCSR<mpf_class> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixCSR<mpf_class> & matrix);
#endif
}
#endif
