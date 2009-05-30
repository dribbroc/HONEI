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

        /// Our row count.
        unsigned long rows;

        /// Our column count.
        unsigned long columns;

        /// Our zero element.
        static const DataType_ zero_element;

        Implementation(unsigned long rows, unsigned long columns) :
            Aj(1),
            Ax(1),
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
            rows(rows),
            columns(columns)
        {
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
        PrivateImplementationPattern<SparseMatrixELL<DataType_>, Shared>(new Implementation<SparseMatrixELL<DataType_> >(src.rows(), src.columns()))
    {
        CONTEXT("When creating SparseMatrixELL from SparseMatrix:");
        unsigned long rows(src.rows());

        unsigned long num_cols_per_row(0);
        for (unsigned long i(0) ; i < rows ; ++i)
        {
            if (src[i].used_elements() > num_cols_per_row)
            {
                num_cols_per_row = src[i].used_elements();
            }
        }
        this->_imp->num_cols_per_row = num_cols_per_row;
        this->_imp->stride = rows;

        DenseVector<unsigned long> Aj(num_cols_per_row * rows, (unsigned long)(0));
        DenseVector<DataType_> Ax(num_cols_per_row * rows, DataType_(0));

        for (unsigned long row(0); row < rows ; ++row)
        {
            DenseVector<DataType_> act_row(src[row]);
            unsigned long i(0);
            unsigned long target(0);
            while(i < act_row.size() && act_row[i] == this->_imp->zero_element)
            {
                ++i;
            }
            for( ; i < act_row.size() ; ++i)
            {
                if (act_row[i] != this->_imp->zero_element)
                {
                    Aj[target + row * num_cols_per_row] = i;
                    Ax[target + row * num_cols_per_row] = act_row[i];
                    target++;
                }
            }
        }

        DenseVector<unsigned long> tAj(num_cols_per_row * rows, (unsigned long)(0));
        DenseVector<DataType_> tAx(num_cols_per_row * rows, DataType_(0));
        unsigned long t(0);
        for (unsigned long j(0) ; j < num_cols_per_row ; ++j)
        {
            for (unsigned long i(0) ; i < tAj.size() ; i+=num_cols_per_row)
            {
                tAj[t] = Aj[j + i];
                tAx[t] = Ax[j + i];
                ++t;
            }
        }
        this->_imp->Aj = tAj;
        this->_imp->Ax = tAx;
    }

    template <typename DataType_>
    SparseMatrixELL<DataType_>::SparseMatrixELL(unsigned long rows, unsigned long columns, DenseVector<unsigned long> & row_indices,
            DenseVector<unsigned long> & column_indices, DenseVector<DataType_> & data) :
        PrivateImplementationPattern<SparseMatrixELL<DataType_>, Shared>(new Implementation<SparseMatrixELL<DataType_> >(rows, columns))
    {
        CONTEXT("When creating SparseMatrixELL from coordinate vectors:");

        SparseMatrix<DataType_> src(rows, columns);
        for (unsigned long i(0) ; i < data.size() ; ++i)
        {
            src(row_indices[i], column_indices[i]) = data[i];
        }

        unsigned long num_cols_per_row(0);
        for (unsigned long i(0) ; i < rows ; ++i)
        {
            if (src[i].used_elements() > num_cols_per_row)
            {
                num_cols_per_row = src[i].used_elements();
            }
        }
        this->_imp->num_cols_per_row = num_cols_per_row;
        this->_imp->stride = rows;

        DenseVector<unsigned long> Aj(num_cols_per_row * rows, (unsigned long)(0));
        DenseVector<DataType_> Ax(num_cols_per_row * rows, DataType_(0));

        for (unsigned long row(0); row < rows ; ++row)
        {
            DenseVector<DataType_> act_row(src[row]);
            unsigned long i(0);
            unsigned long target(0);
            while(i < act_row.size() && act_row[i] == this->_imp->zero_element)
            {
                ++i;
            }
            for( ; i < act_row.size() ; ++i)
            {
                if (act_row[i] != this->_imp->zero_element)
                {
                    Aj[target + row * num_cols_per_row] = i;
                    Ax[target + row * num_cols_per_row] = act_row[i];
                    target++;
                }
            }
        }

        DenseVector<unsigned long> tAj(num_cols_per_row * rows, (unsigned long)(0));
        DenseVector<DataType_> tAx(num_cols_per_row * rows, DataType_(0));
        unsigned long t(0);
        for (unsigned long j(0) ; j < num_cols_per_row ; ++j)
        {
            for (unsigned long i(0) ; i < tAj.size() ; i+=num_cols_per_row)
            {
                tAj[t] = Aj[j + i];
                tAx[t] = Ax[j + i];
                ++t;
            }
        }
        this->_imp->Aj = tAj;
        this->_imp->Ax = tAx;
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
    void SparseMatrixELL<DataType_>::lock(LockMode mode) const
    {
        this->_imp->Aj.lock(mode);
        this->_imp->Ax.lock(mode);
    }

    template <typename DataType_>
            void SparseMatrixELL<DataType_>::unlock(LockMode mode) const
    {
        this->_imp->Aj.unlock(mode);
        this->_imp->Ax.unlock(mode);
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
/*        lhs << "SparseMatrixELL" << std::endl << "{" << std::endl;
        for (unsigned long row(0) ; row < b.size() ; ++row)
        {
            for (unsigned long column(0) ; column < b.size() ; ++column)
            {
                lhs << " " << b(row, column);
            }
            lhs << std::endl;
        }
        lhs << "}" << std::endl;
*/
        lhs << "NumColsPerRow: "<< b.num_cols_per_row() << " Stride: "<< b.stride() << std::endl;
        lhs << "Aj: " << b.Aj();
        lhs << "Ax: " << b.Ax();
        return lhs;
    }
}
#endif
