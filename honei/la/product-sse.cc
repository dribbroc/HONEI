/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
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

#include <honei/la/product.hh>
#include <honei/backends/sse/operations.hh>
#include <honei/la/dense_matrix_tile.hh>

namespace honei
{
    namespace sse
    {
        void rec_dm_product(DenseMatrixTile<float> & r, const DenseMatrixTile<float> & a, const DenseMatrixTile<float> & b)
        {
            /// \todo Use Configuration.
            unsigned long best_cache_size(1000000);
            if (r.rows() * r.columns() + a.rows() * a.columns() + b.rows() * b.columns() < best_cache_size)
            {
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    for (unsigned long j(0) ; j < a.columns() ; ++j)
                    {
                        honei::sse::product_dm(r[i].elements(), b[j].elements(), a(i, j), b[j].size());
                    }
                }
            }
            else
            {
                // case (2.1)
                if (a.rows() >= a.columns() && a.rows() >= b.columns())
                {
                    unsigned long r_lower(r.rows() / 2);
                    unsigned long r_upper(r_lower + r.rows() % 2);
                    DenseMatrixTile<float> r_tile_lower(r, r_lower, r.columns(), 0, 0);
                    DenseMatrixTile<float> r_tile_upper(r, r_upper, r.columns(), r_lower, 0);
                    DenseMatrixTile<float> a_tile_lower(a, r_lower, a.columns(), 0, 0);
                    DenseMatrixTile<float> a_tile_upper(a, r_upper, a.columns(), r_lower, 0);
                    rec_dm_product(r_tile_lower, a_tile_lower, b);
                    rec_dm_product(r_tile_upper, a_tile_upper, b);
                }

                // case (2.3)
                else if(b.columns() >= a.rows() && b.columns() >= a.columns())
                {
                    unsigned long r_left(r.columns() / 2);
                    unsigned long r_right(r_left + r.columns() % 2);
                    DenseMatrixTile<float> r_tile_left(r, r.rows(), r_left, 0, 0);
                    DenseMatrixTile<float> r_tile_right(r, r.rows(), r_right, 0, r_left);
                    DenseMatrixTile<float> b_tile_left(b, b.rows(), r_left, 0, 0);
                    DenseMatrixTile<float> b_tile_right(b, b.rows(), r_right, 0, r_left);
                    rec_dm_product(r_tile_left, a, b_tile_left);
                    rec_dm_product(r_tile_right, a, b_tile_right);
                }

                // case (2.2)
                else
                {
                    unsigned long a_left(a.columns() / 2);
                    unsigned long a_right(a_left + a.columns() % 2);
                    DenseMatrixTile<float> a_tile_left(a, a.rows(), a_left, 0, 0);
                    DenseMatrixTile<float> a_tile_right(a, a.rows(), a_right, 0, a_left);
                    DenseMatrixTile<float> b_tile_upper(b, a_left, b.columns(), 0, 0);
                    DenseMatrixTile<float> b_tile_lower(b, a_right, b.columns(), a_left, 0);
                    // Implicit addition.
                    rec_dm_product(r, a_tile_left, b_tile_upper);
                    rec_dm_product(r, a_tile_right, b_tile_lower);
                }
            }
        }

        void rec_dm_product(DenseMatrixTile<double> & r, const DenseMatrixTile<double> & a, const DenseMatrixTile<double> & b)
        {
            /// \todo Use Configuration.
            unsigned long best_cache_size(500000);
            if (r.rows() * r.columns() + a.rows() * a.columns() + b.rows() * b.columns() < best_cache_size)
            {
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    for (unsigned long j(0) ; j < a.columns() ; ++j)
                    {
                        honei::sse::product_dm(r[i].elements(), b[j].elements(), a(i, j), b[j].size());
                    }
                }
            }
            else
            {
                // case (2.1)
                if (a.rows() >= a.columns() && a.rows() >= b.columns())
                {
                    unsigned long r_lower(r.rows() / 2);
                    unsigned long r_upper(r_lower + r.rows() % 2);
                    DenseMatrixTile<double> r_tile_lower(r, r_lower, r.columns(), 0, 0);
                    DenseMatrixTile<double> r_tile_upper(r, r_upper, r.columns(), r_lower, 0);
                    DenseMatrixTile<double> a_tile_lower(a, r_lower, a.columns(), 0, 0);
                    DenseMatrixTile<double> a_tile_upper(a, r_upper, a.columns(), r_lower, 0);
                    rec_dm_product(r_tile_lower, a_tile_lower, b);
                    rec_dm_product(r_tile_upper, a_tile_upper, b);
                }

                // case (2.3)
                else if(b.columns() >= a.rows() && b.columns() >= a.columns())
                {
                    unsigned long r_left(r.columns() / 2);
                    unsigned long r_right(r_left + r.columns() % 2);
                    DenseMatrixTile<double> r_tile_left(r, r.rows(), r_left, 0, 0);
                    DenseMatrixTile<double> r_tile_right(r, r.rows(), r_right, 0, r_left);
                    DenseMatrixTile<double> b_tile_left(b, b.rows(), r_left, 0, 0);
                    DenseMatrixTile<double> b_tile_right(b, b.rows(), r_right, 0, r_left);
                    rec_dm_product(r_tile_left, a, b_tile_left);
                    rec_dm_product(r_tile_right, a, b_tile_right);
                }

                // case (2.2)
                else
                {
                    unsigned long a_left(a.columns() / 2);
                    unsigned long a_right(a_left + a.columns() % 2);
                    DenseMatrixTile<double> a_tile_left(a, a.rows(), a_left, 0, 0);
                    DenseMatrixTile<double> a_tile_right(a, a.rows(), a_right, 0, a_left);
                    DenseMatrixTile<double> b_tile_upper(b, a_left, b.columns(), 0, 0);
                    DenseMatrixTile<double> b_tile_lower(b, a_right, b.columns(), a_left, 0);
                    // Implicit addition.
                    rec_dm_product(r, a_tile_left, b_tile_upper);
                    rec_dm_product(r, a_tile_right, b_tile_lower);
                }
            }
        }
    }
}

using namespace honei;

DenseVector<float> Product<tags::CPU::SSE>::value(const BandedMatrix<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVectorContinuousBase<float> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows(), float(0));


    unsigned long middle_index(a.rows() - 1);
    unsigned long op_offset;

    for (BandedMatrix<float>::ConstBandIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;
            honei::sse::product_bmdv(result.elements(), band->elements(), b.elements() + op_offset, a.size() - op_offset);
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            honei::sse::product_bmdv(result.elements() + op_offset, band->elements() + op_offset, b.elements(), a.size() - op_offset);
        }
    }

    return result;
}

DenseVector<double> Product<tags::CPU::SSE>::value(const BandedMatrix<double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying BandedMatrix<double> with DenseVectorContinuousBase<double> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows(), double(0));

    unsigned long middle_index(a.rows() - 1);
    unsigned long op_offset;

    for (BandedMatrix<double>::ConstBandIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;
            honei::sse::product_bmdv(result.elements(), band->elements(), b.elements() + op_offset, a.size() - op_offset);
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            honei::sse::product_bmdv(result.elements() + op_offset, band->elements() + op_offset, b.elements(), a.size() - op_offset);
        }
    }

    return result;
}

DenseVector<float> Product<tags::CPU::SSE>::value(const BandedMatrixQ1<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVectorContinuousBase<float> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    honei::sse::product_bmdv_q1(a.band(LL).elements(), a.band(LD).elements(), a.band(LU).elements(),
            a.band(DL).elements(), a.band(DD).elements(), a.band(DU).elements(),
            a.band(UL).elements(), a.band(UD).elements(), a.band(UU).elements(),
            b.elements(), result.elements(), a.size(), a.root());

    return result;
}

DenseVector<double> Product<tags::CPU::SSE>::value(const BandedMatrixQ1<double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying BandedMatrix<double> with DenseVectorContinuousBase<double> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    honei::sse::product_bmdv_q1(a.band(LL).elements(), a.band(LD).elements(), a.band(LU).elements(),
            a.band(DL).elements(), a.band(DD).elements(), a.band(DU).elements(),
            a.band(UL).elements(), a.band(UD).elements(), a.band(UU).elements(),
            b.elements(), result.elements(), a.size(), a.root());

    return result;
}

DenseVector<float> & Product<tags::CPU::SSE>::value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b,
        unsigned long row_start, unsigned long row_end)
{
    CONTEXT("When multiplying SparseMatrixELL<float> with DenseVector<float> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }


    if (row_end == 0)
    {
        fill<tags::CPU::SSE>(result, float(0));
        honei::sse::product_smell_dv(result.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                a.stride(), a.rows(), a.num_cols_per_row());
    }

    else
    {
        honei::sse::product_smell_dv(result.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                a.stride(), a.rows(), a.num_cols_per_row(), row_start, row_end);
    }

    return result;
}

DenseVector<double> & Product<tags::CPU::SSE>::value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b,
        unsigned long row_start, unsigned long row_end)
{
    CONTEXT("When multiplying SparseMatrixELL<double> with DenseVector<double> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    if (row_end == 0)
    {
        fill<tags::CPU::SSE>(result, double(0));
        honei::sse::product_smell_dv(result.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                a.stride(), a.rows(), a.num_cols_per_row());
    }

    else
    {
        honei::sse::product_smell_dv(result.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                a.stride(), a.rows(), a.num_cols_per_row(), row_start, row_end);
    }

    return result;
}

DenseVector<float> Product<tags::CPU::SSE>::value(const DenseMatrix<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> with DenseVectorContinuousBase<float> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    /// \todo Hardcode M*V product.
    for (DenseVector<float>::ElementIterator l(result.begin_elements()), l_end(result.end_elements()) ; l != l_end ; ++l)
    {
        *l = DotProduct<tags::CPU::SSE>::value(b, a[l.index()]);
    }

    return result;
}

DenseVector<double> Product<tags::CPU::SSE>::value(const DenseMatrix<double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> with DenseVectorContinuousBase<double> (SSE):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    /// \todo Hardcode M*V product.
    for (DenseVector<double>::ElementIterator l(result.begin_elements()), l_end(result.end_elements()) ; l != l_end ; ++l)
    {
        *l = DotProduct<tags::CPU::SSE>::value(b, a[l.index()]);
    }

    return result;
}

DenseMatrix<float> Product<tags::CPU::SSE>::value(const DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> with DenseMatrix<float> (SSE):");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<float> result(a.rows(), b.columns(), float(0));

    /*if (a.columns() == a.rows() && b.columns() == 2)
    {
        for (unsigned long i(0) ; i < a.rows() ; ++i)
        {
            honei::sse::product_dm_nx2(result[i].elements(), a[i].elements(), b.elements(), b.rows());
        }
    }
    else
    {*/
        DenseMatrixTile<float> a_tile(a, a.rows(), a.columns(), 0, 0);
        DenseMatrixTile<float> b_tile(b, b.rows(), b.columns(), 0, 0);
        DenseMatrixTile<float> r_tile(result, result.rows(), result.columns(), 0, 0);
        honei::sse::rec_dm_product(r_tile, a_tile, b_tile);
    //}

    return result;
}

DenseMatrix<double> Product<tags::CPU::SSE>::value(const DenseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> with DenseMatrix<double> (SSE):");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<double> result(a.rows(), b.columns(), double(0));

    /*if (a.columns() == a.rows() && b.columns() == 2)
    {
        for (unsigned long i(0) ; i < a.rows() ; ++i)
        {
            honei::sse::product_dm_nx2(result[i].elements(), a[i].elements(), b.elements(), b.rows());
        }
    }
    else
    {*/
        DenseMatrixTile<double> a_tile(a, a.rows(), a.columns(), 0, 0);
        DenseMatrixTile<double> b_tile(b, b.rows(), b.columns(), 0, 0);
        DenseMatrixTile<double> r_tile(result, result.rows(), result.columns(), 0, 0);
        honei::sse::rec_dm_product(r_tile, a_tile, b_tile);
    //}

    return result;
}

DenseMatrix<float> Product<tags::CPU::SSE>::value(const SparseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying SparseMatrix<float> with DenseMatrix<float> (SSE):");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<float> result(a.rows(), b.columns(), float(0));

    for (SparseMatrix<float>::NonZeroConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        honei::sse::product_dm(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}

DenseMatrix<double> Product<tags::CPU::SSE>::value(const SparseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying SparseMatrix<double> with DenseMatrix<double> (SSE):");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<double> result(a.rows(), b.columns(), double(0));

    for (SparseMatrix<double>::NonZeroConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        honei::sse::product_dm(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}

    DenseMatrixTile<float> &
Product<tags::CPU::SSE>::value(DenseMatrixTile<float> & r, const DenseMatrixTile<float> & a, const DenseMatrixTile<float> & b)
{
    CONTEXT("When multiplying DenseMatrixTile<float> with DenseMatrixTile<float> (SSE):");

    if (a.columns() != b.rows() )
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    if (a.rows() != r.rows())
        throw MatrixRowsDoNotMatch(r.rows(), a.rows());

    if (b.columns() != r.columns())
        throw MatrixColumnsDoNotMatch(r.columns(), b.columns());

    honei::sse::rec_dm_product(r, a, b);

    return r;
}

    DenseMatrixTile<double> &
Product<tags::CPU::SSE>::value(DenseMatrixTile<double> & r, const DenseMatrixTile<double> & a, const DenseMatrixTile<double> & b)
{
    CONTEXT("When multiplying DenseMatrixTile<double> with DenseMatrixTile<double> (SSE):");

    if (a.columns() != b.rows() )
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    if (a.rows() != r.rows())
        throw MatrixRowsDoNotMatch(r.rows(), a.rows());

    if (b.columns() != r.columns())
        throw MatrixColumnsDoNotMatch(r.columns(), b.columns());

    honei::sse::rec_dm_product(r, a, b);

    return r;
}
