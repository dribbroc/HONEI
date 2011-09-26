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

#pragma once
#ifndef LIBLA_GUARD_BANDED_MATRIX_QX_HH
#define LIBLA_GUARD_BANDED_MATRIX_QX_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/band_type.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations
    template <typename DataType_> class BandedMatrix;

    /**
     * \brief BandedMatrixQx is a square matrix with nine non-zero bands which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpmatrix
     */
    template <BandType BandType_, typename DataType_> class BandedMatrixQx :
        PrivateImplementationPattern<BandedMatrixQx<BandType_, DataType_>, Shared>
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
            BandedMatrixQx(unsigned long size,
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
            explicit BandedMatrixQx(BandedMatrix<DataType_> & src);

            /**
             * Constructor.
             *
             * \param src The BandedMatrix our matrix will be created from.
             */
            explicit BandedMatrixQx(SparseMatrixELL<DataType_> & src);

            /// Copy-constructor.
            BandedMatrixQx(const BandedMatrixQx<BandType_, DataType_> & other);

            /// Destructor.
            ~BandedMatrixQx();

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
            DenseVector<DataType_> & band(unsigned long index) const;

            /// Returns a band-vector range by index.
            DenseVectorRange<DataType_> band_range(unsigned long index) const;

            /// Request a memory access lock for our data.
            void lock(LockMode mode) const;

            /// Release a memory access lock for our data.
            void unlock(LockMode mode) const;

            /// Returns a copy of the matrix.
            BandedMatrixQx copy() const;

            /// Returns the square root of our size.
            signed long root() const;

            /// Returns number of bands
            unsigned long bands() const;
    };

    /**
     * Equality operator for BandedMatrixQx.
     *
     * Compares if corresponding elements of two banded matrices are equal
     * within machine precision.
     */
    template <BandType BandType_, typename DataType_> bool operator== (const BandedMatrixQx<BandType_, DataType_> & a, const BandedMatrixQx<BandType_, DataType_> & b);

    /**
     * Output operator for BandedMatrixQx.
     *
     * Outputs a banded matrix to an output stream.
     */
    template <BandType BandType_, typename DataType_> std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQx<BandType_, DataType_> & matrix);

    extern template class BandedMatrixQx<Q1Type, float>;

    extern template bool operator== (const BandedMatrixQx<Q1Type, float> & a, const BandedMatrixQx<Q1Type, float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQx<Q1Type, float> & matrix);

    extern template class BandedMatrixQx<Q1Type, double>;

    extern template bool operator== (const BandedMatrixQx<Q1Type, double> & a, const BandedMatrixQx<Q1Type, double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQx<Q1Type, double> & matrix);

}
#endif
