/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DENSE_MATRIX_HH
#define LIBLA_GUARD_DENSE_MATRIX_HH 1

#include <honei/la/const_vector.hh>
#include <honei/la/dense_matrix_tile-fwd.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>
#include <honei/la/element_iterator.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/sparse_matrix-fwd.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    /**
     * DenseMatrix is a matrix with O(column * row) non-zero elements which keeps its data
     * aligned and continuous.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class DenseMatrix :
        public PrivateImplementationPattern<DenseMatrix<DataType_>, Shared>
    {
        public:
            /// \name Our column and row types
            /// \{

            typedef DenseVectorSlice<DataType_> Column;
            typedef DenseVectorSlice<DataType_> ConstColumn;
            typedef DenseVectorRange<DataType_> Row;
            typedef DenseVectorRange<DataType_> ConstRow;

            /// \}

            /// \name Friends of DenseMatrix
            /// \{

            friend class honei::ConstElementIterator<storage::Dense, container::Matrix, DataType_>;
            friend class honei::ElementIterator<storage::Dense, container::Matrix, DataType_>;

            /// \}

            friend class DenseMatrixTile<DataType_>;

            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Dense, container::Matrix, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Dense, container::Matrix, DataType_> ElementIterator;

            /// \name Constructors and destructor
            /// \{

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             */
            DenseMatrix(unsigned long rows, unsigned long columns);

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             * \param value Default value of each of the new matrice's elements.
             */
            DenseMatrix(unsigned long rows, unsigned long columns, DataType_ value);

            /**
             * Constructor.
             *
             * \param other The SparseMatrix to densify.
             */
            DenseMatrix(const SparseMatrix<DataType_> & other);

            /**
             * Constructor
             *
             * Create a submatrix from a given source matrix.
             * \param source The source matrix.
             * \param column_offset The source matrix column offset.
             * \param columns Number of columns of the new matrix.
             * \param row_offset The source matrix row offset.
             * \param rows Number of rows of the new matrix.
             */
            DenseMatrix(const DenseMatrix<DataType_> & source, unsigned long column_offset, unsigned long columns,
                    unsigned long row_offset, unsigned long rows);

            /// Destructor
            ~DenseMatrix();

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const;

            /// Returns iterator pointing to a given element of the matrix.
            virtual ConstElementIterator element_at(unsigned long index) const;

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const;

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements();

            /// Returns iterator pointing to a given element of the matrix.
            virtual ElementIterator element_at(unsigned long index);

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements();

            /// Returns the number of our columns.
            virtual unsigned long columns() const;

            /// Returns the number of our rows.
            virtual unsigned long rows() const;

            /// Returns the number of our elements.
            virtual unsigned long size() const;

            /// Retrieves row vector by index, zero-based, unassignable.
            virtual ConstRow operator[] (unsigned long row) const;

            /// Retrieves row vector by index, zero-based, assignable.
            virtual Row operator[] (unsigned long row);

            /// Retrieves element at (row, column), unassignable.
            virtual const DataType_ & operator() (unsigned long row, unsigned long column) const;

            /// Retrieves element at (row, column), assignable.
            virtual DataType_ & operator() (unsigned long row, unsigned long column);

            /// Retrieves element at (index), unassignable.
            virtual const DataType_ & operator() (unsigned long index) const;

            /// Retrieves element at (index), assignable.
            virtual DataType_ & operator() (unsigned long index);

            /// Retrieves column vector by index, zero-based, unassignable.
            virtual ConstColumn column(unsigned long column) const;

            /// Retrieves column vector by index, zero-based, assignable.
            virtual Column column(unsigned long column);

            /// Returns a pointer to our data array.
            DataType_ * elements() const;

            /// Return the address of our data
            void * address() const;

            /// Return our memory id
            void * memid() const;

            /// Request a memory access lock for our data.
            void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a memory access lock for our data.
            void unlock(LockMode mode) const;

            /// Returns a copy of the matrix.
            DenseMatrix copy() const;
    };

    /**
     * \brief Compare two Matrices for equality.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_> bool operator== (const DenseMatrix<DataType_> & left, const DenseMatrix<DataType_> & right);

    /**
     * \brief Output our DenseMatrix to an ostream.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<DataType_> & m);

    extern template class DenseMatrix<float>;

    extern template bool operator== (const DenseMatrix<float> & a, const DenseMatrix<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<float> & matrix);

    extern template class DenseMatrix<double>;

    extern template bool operator== (const DenseMatrix<double> & a, const DenseMatrix<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<double> & matrix);

    extern template class DenseMatrix<long>;

    extern template bool operator== (const DenseMatrix<long> & a, const DenseMatrix<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<long> & matrix);

    extern template class DenseMatrix<bool>;

    extern template bool operator== (const DenseMatrix<bool> & a, const DenseMatrix<bool> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<bool> & matrix);
}

#endif
