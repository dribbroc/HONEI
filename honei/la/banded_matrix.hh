/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
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

#ifndef LIBLA_GUARD_BANDED_MATRIX_HH
#define LIBLA_GUARD_BANDED_MATRIX_HH 1

#include <honei/la/band_iterator.hh>
#include <honei/la/element_iterator.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix_q1.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations
    template <typename DataType_> class BandedMatrixQ1;

    /**
     * \brief BandedMatrix is a square matrix with O(size) non-zero bands which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class BandedMatrix :
        PrivateImplementationPattern<BandedMatrix<DataType_>, Shared>
    {
        public:
            /// \name Friends of BandedMatrix
            /// \{

            friend struct Implementation<honei::BandIterator<type::Banded, DataType_> >;
            friend struct Implementation<honei::ConstBandIterator<type::Banded, DataType_> >;
            friend struct honei::ConstElementIterator<storage::Banded, container::Matrix, DataType_>;
            friend struct Implementation<honei::ElementIterator<storage::Banded, container::Matrix, DataType_> >;

            /// \}

            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Banded, container::Matrix, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Banded, container::Matrix, DataType_> ElementIterator;

            /// Type of the const iterator over our bands.
            typedef honei::ConstBandIterator<type::Banded, DataType_> ConstBandIterator;

            /// Type of the iterator over our bands.
            typedef honei::BandIterator<type::Banded, DataType_> BandIterator;

            /// \name Basic operations
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             */
            BandedMatrix(unsigned long size);

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             * \param diagonal Diagonal of the new banded matrix.
             */
            BandedMatrix(unsigned long size, const DenseVector<DataType_> & diagonal);

            /**
             * Constructor.
             *
             * \param source The q1 matrix our matrix will be created from.
             */
            explicit BandedMatrix(const BandedMatrixQ1<DataType_> & src);

            /// Destructor.
            ~BandedMatrix();

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the matrix.
            ElementIterator end_elements();

            /// Returns const iterator pointing to the first element of the matrix.
            ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the matrix.
            ConstElementIterator end_elements() const;

            /// Returns iterator pointing to the first band of the matrix.
            BandIterator begin_bands();

            /// Returns iterator pointing to a given band of the matrix.
            BandIterator band_at(unsigned long index);

            /// Returns iterator pointing behind the last band of the matrix.
            BandIterator end_bands();

            /// Returns iterator pointing to the first band of the matrix.
            ConstBandIterator begin_bands() const;

            /// Returns iterator pointing to a given band of the matrix.
            ConstBandIterator band_at(unsigned long index) const;

            /// Returns iterator pointing behind the last band of the matrix.
            ConstBandIterator end_bands() const;

            /// Returns iterator pointing to the first inserted non zero band of the matrix.
            BandIterator begin_non_zero_bands();

            /// Returns iterator pointing behind the last inserted band of the matrix.
            BandIterator end_non_zero_bands();

            /// Returns iterator pointing to the first inserted non zero band of the matrix.
            ConstBandIterator begin_non_zero_bands() const;

            /// Returns iterator pointing behind the last inserted band of the matrix.
            ConstBandIterator end_non_zero_bands() const;

            /// Returns the number of our columns.
            unsigned long columns() const;

            /// Returns the number of our rows.
            unsigned long rows() const;

            /// Returns our size, equal to rows and columns.
            unsigned long size() const;

            /// Returns the number of our non zero bands.
            unsigned long non_zero_bands() const;

            /// Inserts a new Band in the matrix.
            void insert_band(signed long index, const DenseVector<DataType_> & vector);

            /// Returns a band-vector by index.
            DenseVector<DataType_> & band(signed long index) const;

            /// Returns a band-vector by unsigned index.
            DenseVector<DataType_> & band_unsigned(unsigned long index);

            /// Request a memory access lock for our data.
            void lock(LockMode mode) const;

            /// Release a memory access lock for our data.
            void unlock(LockMode mode) const;

            /// Returns a copy of the matrix.
            BandedMatrix copy() const;
    };

    /**
     * Equality operator for BandedMatrix.
     *
     * Compares if corresponding elements of two banded matrices are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const BandedMatrix<DataType_> & a, const BandedMatrix<DataType_> & b);

    /**
     * Output operator for BandedMatrix.
     *
     * Outputs a banded matrix to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const BandedMatrix<DataType_> & matrix);

    extern template class BandedMatrix<float>;

    extern template bool operator== (const BandedMatrix<float> & a, const BandedMatrix<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrix<float> & matrix);

    extern template class BandedMatrix<double>;

    extern template bool operator== (const BandedMatrix<double> & a, const BandedMatrix<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrix<double> & matrix);
}
#endif
