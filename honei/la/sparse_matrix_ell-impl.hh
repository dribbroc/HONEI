/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LIBLA_GUARD_SPARSE_MATRIX_ELL_IMPL_HH
#define LIBLA_GUARD_SPARSE_MATRIX_ELL_IMPL_HH 1

#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/assertion.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <cmath>

namespace honei
{
    template <typename DataType_> struct Implementation<SparseMatrixELL<DataType_> >
    {
        unsigned long stride;
        unsigned long num_cols_per_row;

        DenseVector<unsigned long> Aj;//column indices stored in a (cols_per_row x stride) matrix
        DenseVector<DataType_> Ax;//nonzero values stored in a (cols_per_row x stride) matrix
        DenseVector<unsigned long> Arl;//length of every single row

        /// Our row count.
        unsigned long rows;

        /// Our column count.
        unsigned long columns;

        /// Our zero element.
        static const DataType_ zero_element;

        Implementation(unsigned long rows, unsigned long columns) :
            Aj(1),
            Ax(1),
            Arl(1),
            rows(rows),
            columns(columns)
        {
        }

        Implementation(unsigned long rows, unsigned long columns, unsigned long stride, unsigned long num_cols_per_row,
                const DenseVector<unsigned long> & Aj, const DenseVector<DataType_> & Ax) :
            stride(stride),
            num_cols_per_row(num_cols_per_row),
            Aj(Aj),
            Ax(Ax),
            Arl(1),
            rows(rows),
            columns(columns)
        {
            Arl = row_length();
        }

        Implementation(SparseMatrix<DataType_> & src) :
            Aj(1),
            Ax(1),
            Arl(1),
            rows(1),
            columns(1)
        {
            rows = src.rows();
            columns = src.columns();

            num_cols_per_row = 1;
            for (unsigned long i(0) ; i < rows ; ++i)
            {
                if (src[i].used_elements() > num_cols_per_row)
                {
                    num_cols_per_row = src[i].used_elements();
                }
            }
            /// \todo remove hardcoded numbers
            unsigned long alignment(32);
            stride = alignment * ((rows + alignment - 1)/ alignment);

            DenseVector<unsigned long> pAj(num_cols_per_row * stride, (unsigned long)(0));
            DenseVector<DataType_> pAx(num_cols_per_row * stride, DataType_(0));

            for (unsigned long row(0); row < rows ; ++row)
            {
                unsigned long target(0);
                for (typename SparseVector<DataType_>::NonZeroElementIterator i(src[row].begin_non_zero_elements()) ;
                        i < src[row].end_non_zero_elements() ; ++i)
                {
                    if(*i != DataType_(0))
                    {
                        pAj[row + target * stride] = i.index();
                        pAx[row + target * stride] = *i;
                        target++;
                    }
                }
            }

            Aj = pAj;
            Ax = pAx;
            Arl = row_length();
        }

        private:
        DenseVector<unsigned long> row_length()
        {
            DenseVector<unsigned long> rl(rows, 0);
            for (unsigned long row(0) ; row < rows ; ++row)
            {
                unsigned long count(0);
                for (unsigned long i(row) ; i < Ax.size() ; i+=stride)
                {
                    if (Ax[i] == DataType_(0))
                        break;
                    ++count;
                }
                rl[row] = count;
            }
            return rl;
        }
    };


    template <typename DataType_>
    SparseMatrixELL<DataType_>::SparseMatrixELL(unsigned long rows, unsigned columns, unsigned long stride,
            unsigned long num_cols_per_row,
            const DenseVector<unsigned long> & Aj,
            const DenseVector<DataType_> & Ax) :
        PrivateImplementationPattern<SparseMatrixELL<DataType_>, Shared>(new Implementation<SparseMatrixELL<DataType_> >(rows, columns, stride, num_cols_per_row, Aj, Ax))
    {
        CONTEXT("When creating SparseMatrixELL:");
    }

    template <typename DataType_>
    SparseMatrixELL<DataType_>::SparseMatrixELL(SparseMatrix<DataType_> & src) :
        PrivateImplementationPattern<SparseMatrixELL<DataType_>, Shared>(new Implementation<SparseMatrixELL<DataType_> >(src))
    {
        CONTEXT("When creating SparseMatrixELL from SparseMatrix:");
    }

    template <typename DataType_>
    SparseMatrixELL<DataType_>::SparseMatrixELL(const SparseMatrixELL<DataType_> & other) :
        PrivateImplementationPattern<SparseMatrixELL<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    SparseMatrixELL<DataType_>::~SparseMatrixELL()
    {
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::columns() const
    {
        return this->_imp->columns;
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::rows() const
    {
        return this->_imp->rows;
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::size() const
    {
        return this->_imp->columns * this->_imp->rows;
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::used_elements() const
    {
        unsigned long ue(0);
        for(unsigned long i(0) ; i < this->_imp->Arl.size() ; ++i)
            ue += this->_imp->Arl.elements()[i];

        return ue;
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::stride() const
    {
        return this->_imp->stride;
    }

    template <typename DataType_>
    unsigned long
    SparseMatrixELL<DataType_>::num_cols_per_row() const
    {
        return this->_imp->num_cols_per_row;
    }

    template <typename DataType_>
    DenseVector<unsigned long> &
    SparseMatrixELL<DataType_>::Aj() const
    {
        return this->_imp->Aj;
    }

    template <typename DataType_>
    DenseVector<DataType_> &
    SparseMatrixELL<DataType_>::Ax() const
    {
        return this->_imp->Ax;
    }

    template <typename DataType_>
    DenseVector<unsigned long> &
    SparseMatrixELL<DataType_>::Arl() const
    {
        return this->_imp->Arl;
    }

    template <typename DataType_>
    const DataType_ SparseMatrixELL<DataType_>::operator() (unsigned long row, unsigned long column) const
    {
        unsigned long max(this->Arl()[row]);
        for (unsigned long i(row), j(0) ; j < max && this->Aj()[i] <= column ; i += this->stride(), ++j)
        {
            if (this->Aj()[i] == column)
                return this->Ax()[i];
        }
        return this->_imp->zero_element;
    }

    template <typename DataType_>
    void SparseMatrixELL<DataType_>::lock(LockMode mode) const
    {
        this->_imp->Aj.lock(mode);
        this->_imp->Ax.lock(mode);
        this->_imp->Arl.lock(mode);
    }

    template <typename DataType_>
            void SparseMatrixELL<DataType_>::unlock(LockMode mode) const
    {
        this->_imp->Aj.unlock(mode);
        this->_imp->Ax.unlock(mode);
        this->_imp->Arl.unlock(mode);
    }

    template <typename DataType_>
    SparseMatrixELL<DataType_>
    SparseMatrixELL<DataType_>::copy() const
    {
        CONTEXT("When creating copy() of a SparseMatrixELL:");
        SparseMatrixELL result(this->_imp->rows,
                this->_imp->columns,
                this->_imp->stride,
                this->_imp->num_cols_per_row,
                this->_imp->Aj.copy(),
                this->_imp->Ax.copy());

        return result;
    }

    template <typename DataType_>
    bool
    operator== (const SparseMatrixELL<DataType_> & a, const SparseMatrixELL<DataType_> & b)
    {
        if (a.columns() != b.columns())
        {
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
        }

        if (a.rows() != b.rows())
        {
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        }

        bool result(true);

        result &= (a.Aj() == b.Aj());
        result &= (a.Ax() == b.Ax());

        return result;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrixELL<DataType_> & b)
    {
        lhs << "SparseMatrixELL" << std::endl << "[" << std::endl;
        for (unsigned long row(0) ; row < b.rows() ; ++row)
        {
            lhs << "[";
            for (unsigned long column(0) ; column < b.columns() ; ++column)
            {
                lhs << " " << b(row, column);
            }
            lhs << "]";
            lhs << std::endl;
        }
        lhs << "]" << std::endl;

        lhs << "NumColsPerRow: "<< b.num_cols_per_row() << " Stride: "<< b.stride() << std::endl;
        lhs << "Aj: " << b.Aj();
        lhs << "Ax: " << b.Ax();
        lhs << "Arl: " << b.Arl();
        return lhs;
    }
}
#endif
