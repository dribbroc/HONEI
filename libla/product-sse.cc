
/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <libla/product.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            extern void scaled_sum(float * x, const float * y, float b, unsigned long size);
            extern void scaled_sum(double * x, const double * y, double b, unsigned long size);
        }
    }
}

using namespace honei;

DenseVector<float> Product<tags::CPU::SSE>::value(const BandedMatrix<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVector<float> with SSE:");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows(), float(0));

    __m128 m1, m2, m3, m4, m5, m6;

    unsigned long middle_index(a.rows() - 1);
    unsigned long quad_end, end, quad_start, start, op_offset;

    // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
    for (BandedMatrix<float>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
            vi != vi_end ; ++vi)
    {
        if (! vi.exists())
            continue;
        op_offset = vi.index() - middle_index;
        end = vi->size() - op_offset; //Calculation of the element-index to stop in iteration!
        quad_end = end - (end % 8);
        if (end < 32) quad_end = 0;

        float * band_e = vi->elements();
        float * b_e = b.elements();
        float * r_e = result.elements();

        for (unsigned long index = 0 ; index < quad_end ; index += 8)
        {
            m2 = _mm_loadu_ps(b_e + index + op_offset);
            m5 = _mm_loadu_ps(b_e + index + op_offset + 4);
            m1 = _mm_load_ps(band_e + index);
            m4 = _mm_load_ps(band_e + index + 4);
            m3 = _mm_load_ps(r_e + index);
            m6 = _mm_load_ps(r_e + index + 4);

            m1 = _mm_mul_ps(m1, m2);
            m4 = _mm_mul_ps(m4, m5);
            m1 = _mm_add_ps(m1, m3);
            m4 = _mm_add_ps(m4, m6);

            _mm_store_ps(r_e + index, m1);
            _mm_store_ps(r_e + index + 4, m4);
        }

        for (unsigned long index = quad_end ; index < end ; index++) 
        {
            r_e[index] += band_e[index] * b_e[index + op_offset];
        }
    }

    // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
    for (BandedMatrix<float>::ConstVectorIterator vi(a.begin_bands()), vi_end(a.band_at(middle_index)) ;
            vi != vi_end ; ++vi)
    {
        if (! vi.exists())
            continue;
        op_offset = middle_index - vi.index();
        start = op_offset; //Calculation of the element-index to start in iteration!
        quad_start = start + (8 - (start % 8));
        end = a.size();
        quad_end = end - (end % 8);
        if ( start + 32 > end)
        {
            quad_end = start;
            quad_start = start;
        }
        float * band_e = vi->elements();
        float * b_e = b.elements();
        float * r_e = result.elements();

        for (unsigned long index = start ; index < quad_start ; index++)
        {
            r_e[index] += band_e[index] * b_e[index - op_offset];
        }
        for (unsigned long index = quad_start ; index < quad_end ; index += 8)
        {
            m2 = _mm_loadu_ps(b_e + index - op_offset);
            m5 = _mm_loadu_ps(b_e + index - op_offset + 4);
            m1 = _mm_load_ps(band_e + index);
            m4 = _mm_load_ps(band_e + index + 4);
            m3 = _mm_load_ps(r_e + index);
            m6 = _mm_load_ps(r_e + index + 4);

            m1 = _mm_mul_ps(m1, m2);
            m4 = _mm_mul_ps(m4, m5);
            m1 = _mm_add_ps(m1, m3);
            m4 = _mm_add_ps(m4, m6);

            _mm_store_ps(r_e + index, m1);
            _mm_store_ps(r_e + index + 4, m4);
        }

        for (unsigned long index = quad_end ; index < end ; index++)
        {
            r_e[index] += band_e[index] * b_e[index - op_offset];
        }
    }
    return result;
}

DenseVector<double> Product<tags::CPU::SSE>::value(const BandedMatrix<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying BandedMatrix<double> with DenseVector<double> with SSE:");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows(), double(0));

    __m128d m1, m2, m3, m4, m5, m6;

    unsigned long middle_index(a.rows() - 1);
    unsigned long quad_end, end, quad_start, start, op_offset;

    // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
    for (BandedMatrix<double>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
            vi != vi_end ; ++vi)
    {
        if (! vi.exists())
            continue;
        op_offset = vi.index() - middle_index;
        end = vi->size() - op_offset; //Calculation of the element-index to stop in iteration!
        quad_end = end - (end % 4);
        if (end < 12) quad_end = 0;

        double * band_e = vi->elements();
        double * b_e = b.elements();
        double * r_e = result.elements();

        for (unsigned long index = 0 ; index < quad_end ; index += 4)
        {
            m2 = _mm_loadu_pd(b_e + index + op_offset);
            m5 = _mm_loadu_pd(b_e + index + op_offset + 2);
            m1 = _mm_load_pd(band_e + index);
            m4 = _mm_load_pd(band_e + index + 2);
            m3 = _mm_load_pd(r_e + index);
            m6 = _mm_load_pd(r_e + index + 2);

            m1 = _mm_mul_pd(m1, m2);
            m4 = _mm_mul_pd(m4, m5);
            m1 = _mm_add_pd(m1, m3);
            m4 = _mm_add_pd(m4, m6);

            _mm_store_pd(r_e + index, m1);
            _mm_store_pd(r_e + index + 2, m4);
        }

        for (unsigned long index = quad_end ; index < end ; index++) 
        {
            r_e[index] += band_e[index] * b_e[index + op_offset];
        }
    }

    // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
    for (BandedMatrix<double>::ConstVectorIterator vi(a.begin_bands()), vi_end(a.band_at(middle_index)) ;
            vi != vi_end ; ++vi)
    {
        if (! vi.exists())
            continue;
        op_offset = middle_index - vi.index();
        start = op_offset; //Calculation of the element-index to start in iteration!
        quad_start = start + (4 - (start % 4));
        end = a.size();
        quad_end = end - (end % 4);
        if ( start + 12 > end)
        {
            quad_end = start;
            quad_start = start;
        }
        double * band_e = vi->elements();
        double * b_e = b.elements();
        double * r_e = result.elements();

        for (unsigned long index = start ; index < quad_start ; index++)
        {
            r_e[index] += band_e[index] * b_e[index - op_offset];
        }
        for (unsigned long index = quad_start ; index < quad_end ; index += 4)
        {
            m2 = _mm_loadu_pd(b_e + index - op_offset);
            m5 = _mm_loadu_pd(b_e + index - op_offset + 2);
            m1 = _mm_load_pd(band_e + index);
            m4 = _mm_load_pd(band_e + index + 2);
            m3 = _mm_load_pd(r_e + index);
            m6 = _mm_load_pd(r_e + index + 2);

            m1 = _mm_mul_pd(m1, m2);
            m4 = _mm_mul_pd(m4, m5);
            m1 = _mm_add_pd(m1, m3);
            m4 = _mm_add_pd(m4, m6);

            _mm_store_pd(r_e + index, m1);
            _mm_store_pd(r_e + index + 2, m4);
        }

        for (unsigned long index = quad_end ; index < end ; index++)
        {
            r_e[index] += band_e[index] * b_e[index - op_offset];
        }
    }
    return result;
}

DenseVector<float> Product<tags::CPU::SSE>::value(const DenseMatrix<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> with DenseVector<float> with SSE:");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    /// \todo Hardcode M*V product.
    for (Vector<float>::ElementIterator l(result.begin_elements()), l_end(result.end_elements()) ; l != l_end ; ++l)
    {
        *l = DotProduct<tags::CPU::SSE>::value(b, a[l.index()]);
    }

    return result;
}

DenseVector<double> Product<tags::CPU::SSE>::value(const DenseMatrix<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> with DenseVector<double> with SSE:");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows());

    /// \todo Hardcode M*V product.
    for (Vector<double>::ElementIterator l(result.begin_elements()), l_end(result.end_elements()) ; l != l_end ; ++l)
    {
        *l = DotProduct<tags::CPU::SSE>::value(b, a[l.index()]);
    }

    return result;
}

DenseMatrix<float> Product<tags::CPU::SSE>::value(const DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> with DenseMatrix<float> with SSE:");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<float> result(a.rows(), b.columns(), float(0));

    for (DenseMatrix<float>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()) ;
            i != i_end ; ++i)
    {
        honei::intern::sse::scaled_sum(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}

DenseMatrix<double> Product<tags::CPU::SSE>::value(const DenseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> with DenseMatrix<double> with SSE:");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<double> result(a.rows(), b.columns(), double(0));

    for (DenseMatrix<double>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()) ;
            i != i_end ; ++i)
    {
        honei::intern::sse::scaled_sum(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}

DenseMatrix<float> Product<tags::CPU::SSE>::value(const SparseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying SparseMatrix<float> with DenseMatrix<float> with SSE:");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<float> result(a.rows(), b.columns(), float(0));

    for (SparseMatrix<float>::ConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        honei::intern::sse::scaled_sum(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}

DenseMatrix<double> Product<tags::CPU::SSE>::value(const SparseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying SparseMatrix<double> with DenseMatrix<double> with SSE:");

    if (a.columns() != b.rows())
        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

    DenseMatrix<double> result(a.rows(), b.columns(), double(0));

    for (SparseMatrix<double>::ConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        honei::intern::sse::scaled_sum(result[i.row()].elements(), b[i.column()].elements(), *i, b[i.column()].size());
    }

    return result;
}
