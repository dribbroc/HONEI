/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/libmath/jacobi_kernel.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

DenseVector<float> JacobiKernel<tags::CPU::SSE>::value(DenseVector<float> & b, DenseVector<float> & x, DenseVector<float> & d, BandedMatrix<float> & a)
{

    CONTEXT("When performing singlestep Jacobi method with SSE.");
    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows(), float(0));

    __m128 m1, m2, m3, m4, m5, m6, m7, m8; //Two in addition for the diag_inv and rhs vectors!!!

    unsigned long middle_index(a.rows() - 1);
    unsigned long quad_end, end, quad_start, start, op_offset;
    DenseVector<float> temp(result.copy());
    float * r_e_c = temp.elements();
    float * x_e = x.elements();
    float * r_e = result.elements();
    float * d_e = d.elements();

    for (BandedMatrix<float>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {

        float * band_e = band->elements();
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {

            op_offset = band.index() - middle_index;

            end = a.size() - op_offset; // Calculation of the element-index to stop in iteration!
            quad_end = end - (end % 4);

            if (end < 32)
                quad_end = 0;

            for (unsigned long index(0) ; index < quad_end ; index += 4)
            {
                m2 = _mm_loadu_ps(x_e + index + op_offset);
                m1 = _mm_load_ps(band_e + index);
                m3 = _mm_load_ps(r_e_c + index);
                m4 = _mm_load_ps(d_e + index);
                m6 = _mm_load_ps(r_e + index);

                m1 = _mm_mul_ps(m1, m2);
                m5 = _mm_add_ps(m1, m3);

                _mm_store_ps(r_e_c + index, m5);

                m1 = _mm_mul_ps(m1, m4);
                m1 = _mm_add_ps(m1, m6);

                _mm_store_ps(r_e + index, m1);
            }

            for (unsigned long index(quad_end) ; index < end ; index++) 
            {
                r_e_c[index] += (band_e[index] * x_e[index + op_offset]);
                r_e[index] += (band_e[index] * x_e[index + op_offset] * d_e[index]);
            }
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();

            start = op_offset; // Calculation of the element-index to start in iteration!
            quad_start = start + (4 - (start % 4));
            end = a.size();
            quad_end = end - (end % 4);

            if ( start + 32 > end)
            {
                quad_end = start;
                quad_start = start;
            }


            for (unsigned long index(start) ; index < quad_start ; index++)
            {
                r_e_c[index] += (band_e[index] * x_e[index - op_offset]);
                r_e[index] += (band_e[index] * x_e[index - op_offset] * d_e[index]);
            }

            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                m2 = _mm_loadu_ps(x_e + index - op_offset);
                m1 = _mm_load_ps(band_e + index);
                m3 = _mm_load_ps(r_e_c + index);
                m4 = _mm_load_ps(d_e + index);
                m6 = _mm_load_ps(r_e + index);

                m1 = _mm_mul_ps(m1, m2);
                m5 = _mm_add_ps(m1, m3);

                _mm_store_ps(r_e_c + index, m5);

                m1 = _mm_mul_ps(m1, m4);
                m1 = _mm_add_ps(m1, m6);

                _mm_store_ps(r_e + index, m1);

            }
            for (unsigned long index(quad_end) ; index < end ; index++)
            {
                r_e_c[index] += (band_e[index] * x_e[index - op_offset]);
                r_e[index] +=  (band_e[index] * x_e[index - op_offset] * d_e[index]);
            }
        }

    }
    //treat the d_e *.elem b_e product as if it was an additional band:
    float * b_e = b.elements();
    end = b.size(); // Calculation of the element-index to stop in iteration!
    quad_end = end - (end % 4);

    if (end < 32)
        quad_end = 0;

    for (unsigned long index(0) ; index < quad_end ; index += 4)
    {
        m1 = _mm_load_ps(d_e + index);
        m2 = _mm_load_ps(r_e + index);
        m3 = _mm_load_ps(b_e + index);
        m1 = _mm_mul_ps(m1, m3);
        m1 = _mm_add_ps(m1, m2);

        _mm_store_ps(r_e + index, m1);
    }

    for (unsigned long index(quad_end) ; index < end ; index++) 
    {
        r_e[index] += (d_e[index] * b_e[index]);
    }

    x = result.copy();
    return result;
}

DenseVector<double> JacobiKernel<tags::CPU::SSE>::value(DenseVector<double> & b, DenseVector<double> & x, DenseVector<double> & d, BandedMatrix<double> & a)
{

    CONTEXT("When performing singlestep Jacobi method with SSE.");
    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<double> result(a.rows(), double(0));

    __m128d m1, m2, m3, m4, m5, m6, m7, m8; //Two in addition for the diag_inv and rhs vectors!!!

    unsigned long middle_index(a.rows() - 1);
    unsigned long quad_end, end, quad_start, start, op_offset;
    DenseVector<double> temp(result.copy());
    double * r_e_c = temp.elements();
    double * x_e = x.elements();
    double * r_e = result.elements();
    double * d_e = d.elements();
    double * b_e = b.elements();

    for (BandedMatrix<double>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;

            end = a.size() - op_offset; // Calculation of the element-index to stop in iteration!
            quad_end = end - (end % 2);

            if (end < 32)
                quad_end = 0;

            double * band_e = band->elements();

            for (unsigned long index(0) ; index < quad_end ; index += 2)
            {
                m2 = _mm_loadu_pd(x_e + index + op_offset);
                m1 = _mm_load_pd(band_e + index);
                m3 = _mm_load_pd(r_e_c + index);
                m4 = _mm_load_pd(d_e + index);
                m6 = _mm_load_pd(r_e + index);

                m1 = _mm_mul_pd(m1, m2);
                m5 = _mm_add_pd(m1, m3);

                _mm_store_pd(r_e_c + index, m5);

                m1 = _mm_mul_pd(m1, m4);
                m1 = _mm_add_pd(m1, m6);

                _mm_store_pd(r_e + index, m1);
            }

            for (unsigned long index(quad_end) ; index < end ; index++) 
            {
                r_e_c[index] += (band_e[index] * x_e[index + op_offset]);
                r_e[index] += (band_e[index] * x_e[index + op_offset] * d_e[index]);
            }
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();

            start = op_offset; // Calculation of the element-index to start in iteration!
            quad_start = start + (2 - (start % 2));
            end = a.size();
            quad_end = end - (end % 2);

            if ( start + 32 > end)
            {
                quad_end = start;
                quad_start = start;
            }

            double * band_e = band->elements();

            for (unsigned long index(start) ; index < quad_start ; index++)
            {
                r_e_c[index] += (band_e[index] * x_e[index - op_offset]);
                r_e[index] += (band_e[index] * x_e[index - op_offset] * d_e[index]);
            }

            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                m2 = _mm_loadu_pd(x_e + index - op_offset);
                m1 = _mm_load_pd(band_e + index);
                m3 = _mm_load_pd(r_e_c + index);
                m4 = _mm_load_pd(d_e + index);
                m6 = _mm_load_pd(r_e + index);

                m1 = _mm_mul_pd(m1, m2);
                m5 = _mm_add_pd(m1, m3);

                _mm_store_pd(r_e_c + index, m5);

                m1 = _mm_mul_pd(m1, m4);
                m1 = _mm_add_pd(m1, m6);

                _mm_store_pd(r_e + index, m1);

            }
            for (unsigned long index(quad_end) ; index < end ; index++)
            {
                r_e_c[index] += (band_e[index] * x_e[index - op_offset]);
                r_e[index] +=  (band_e[index] * x_e[index - op_offset] * d_e[index]);
            }
        }

    }
    //treat the d_e *.elem b_e product as if it was an additional band:
    end = b.size(); // Calculation of the element-index to stop in iteration!
    quad_end = end - (end % 2);

    if (end < 32)
        quad_end = 0;

    for (unsigned long index(0) ; index < quad_end ; index += 2)
    {
        m1 = _mm_load_pd(d_e + index);
        m2 = _mm_load_pd(r_e + index);
        m3 = _mm_load_pd(b_e + index);
        m1 = _mm_mul_pd(m1, m3);
        m1 = _mm_add_pd(m1, m2);

        _mm_store_pd(r_e + index, m1);
    }

    for (unsigned long index(quad_end) ; index < end ; index++) 
    {
        r_e[index] += (d_e[index] * b_e[index]);
    }

    x = result.copy();
    return result;
}
