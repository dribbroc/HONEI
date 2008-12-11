/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_BANDED_MATRIX_Q1_HH
#define LIBLA_GUARD_BANDED_MATRIX_Q1_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations
    template <typename DataType_> class BandedMatrix;

    enum Q1BandIndex
    {
        LL = 0,
        LD,
        LU,
        DL,
        DD,
        DU,
        UL,
        UD,
        UU
    };

    /**
     * \brief BandedMatrixQ1 is a square matrix with nine non-zero bands which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class BandedMatrixQ1 :
        PrivateImplementationPattern<BandedMatrixQ1<DataType_>, Shared>
    {
        public:
            /// \name Basic operations
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new banded matrix.
             * \param diagonal Diagonal of the new banded matrix.
             */
            BandedMatrixQ1(unsigned long size,
                    const DenseVector<DataType_> & ll,
                    const DenseVector<DataType_> & ld,
                    const DenseVector<DataType_> & lu,
                    const DenseVector<DataType_> & dl,
                    const DenseVector<DataType_> & dd,
                    const DenseVector<DataType_> & du,
                    const DenseVector<DataType_> & ul,
                    const DenseVector<DataType_> & ud,
                    const DenseVector<DataType_> & uu);

            /**
             * Constructor.
             *
             * \param src The BandedMatrix our matrix will be created from.
             */
            BandedMatrixQ1(BandedMatrix<DataType_> & src);

            /// Copy-constructor.
            BandedMatrixQ1(const BandedMatrixQ1<DataType_> & other);

            /// Destructor.
            ~BandedMatrixQ1();

            /// \}

            /// Returns the number of our columns.
            unsigned long columns() const;

            /// Returns the number of our rows.
            unsigned long rows() const;

            /// Returns our size, equal to rows and columns.
            unsigned long size() const;

            /// Retrieves element at (row, column), unassignable.
            const DataType_ & operator() (unsigned long row, unsigned long column) const;

            /// Returns a full band-vector by index.
            DenseVector<DataType_> & band(Q1BandIndex index) const;

            /// Returns a band-vector range by index.
            DenseVectorRange<DataType_> band_range(Q1BandIndex index) const;

            /// Request a memory access lock for our data.
            void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a memory access lock for our data.
            void unlock(LockMode mode) const;

            /// Returns a copy of the matrix.
            BandedMatrixQ1 copy() const;

            /// Returns the square root of our size.
            signed long root() const;
    };

    /**
     * Equality operator for BandedMatrixQ1.
     *
     * Compares if corresponding elements of two banded matrices are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const BandedMatrixQ1<DataType_> & a, const BandedMatrixQ1<DataType_> & b);

    /**
     * Output operator for BandedMatrixQ1.
     *
     * Outputs a banded matrix to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQ1<DataType_> & matrix);

    extern template class BandedMatrixQ1<float>;

    extern template bool operator== (const BandedMatrixQ1<float> & a, const BandedMatrixQ1<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQ1<float> & matrix);

    extern template class BandedMatrixQ1<double>;

    extern template bool operator== (const BandedMatrixQ1<double> & a, const BandedMatrixQ1<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQ1<double> & matrix);

}
#endif
