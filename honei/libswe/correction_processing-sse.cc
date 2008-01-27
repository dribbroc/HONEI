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

#include <honei/libswe/correction_processing.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

void CorrectionProcessing<boundaries::REFLECT, tags::CPU::SSE>
                        ::value(DenseVector<float> & predictedu,
                                DenseVector<float> & predictedv,
                                DenseVector<float> & predictedw,
                                DenseVector<float> & u,
                                DenseVector<float> & v,
                                DenseVector<float> & w,
                                unsigned long d_width,
                                unsigned long d_height,
                                DenseMatrix<float> & height)
{
    float * u_e = u.elements();
    float * v_e = v.elements();
    float * w_e = w.elements();

    float * u_p_e = predictedu.elements();
    float * v_p_e = predictedv.elements();
    float * w_p_e = predictedw.elements();

    float * h_e = height.elements();

    unsigned long sse_steps((u.size() - (u.size() % 4)) / 4);
    unsigned long sse_limit(sse_steps * 4);
    ///SSE loop:
    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float __attribute__((aligned(16))) one_half(0.5f);
    m1 = _mm_set_ps1(one_half);


    float __attribute__((aligned(16))) v_c_e[v.size()];
    float __attribute__((aligned(16))) u_c_e[u.size()];
    float __attribute__((aligned(16))) w_c_e[w.size()];

    unsigned long c_count(0);
    for(unsigned long i(0); i < u.size(); ++i)
    {
        if(i < 6*(d_width+4)+6 || i >3 * ((d_height + 2) * (d_width + 4)))
        {

            v_c_e[i] = v_e[i];
            u_c_e[i] = u_e[i];
            w_c_e[i] = w_e[i];

            ++c_count;
        }
        else if(c_count == d_width)
        {
            for(unsigned j(0); j < 12; ++j)
            {
                u_c_e[i + j] = u_e[i + j];
                v_c_e[i + j] = v_e[i + j];
                w_c_e[i + j] = w_e[i + j];
            }
            c_count =0;
            i += 11;
        }

    }

    for(unsigned long i(0); i < sse_limit; i += 4)
    {
        m2 = _mm_load_ps(u_e + i);
        m3 = _mm_load_ps(v_e + i);
        m4 = _mm_load_ps(w_e + i);
        m5 = _mm_load_ps(u_p_e + i);
        m6 = _mm_load_ps(v_p_e + i);
        m7 = _mm_load_ps(w_p_e + i);

        m8 = _mm_add_ps(m2, m5);
        m2 = _mm_mul_ps(m8, m1);

        m8 = _mm_add_ps(m3, m6);
        m3 = _mm_mul_ps(m8, m1);

        m8 = _mm_add_ps(m4, m7);
        m4 = _mm_mul_ps(m8, m1);

        _mm_store_ps(u_e + i, m2);
        _mm_store_ps(v_e + i, m3);
        _mm_store_ps(w_e + i, m4);
    }

    ///Non SSE computation:
    for(unsigned long i(sse_limit); i < u.size(); ++i)
    {
        u_e[i] = float(0.5) * (u_e[i] + u_p_e[i]);
        v_e[i] = float(0.5) * (v_e[i] + v_p_e[i]);
        w_e[i] = float(0.5) * (w_e[i] + w_p_e[i]);
    }

    unsigned long col_count(0);
    for(unsigned long i(0); i < u.size(); ++i)
    {
        if(i < 6*(d_width+4)+6 || i >3 * ((d_height + 2) * (d_width + 4)))
        {
            u_e[i] = u_c_e[i];
            v_e[i] = v_c_e[i];
            w_e[i] = w_c_e[i];

            ++col_count;
        }
        else if(col_count == d_width)
        {
            for(unsigned j(0); j < 12; ++j)
            {
                u_e[i + j] = u_c_e[i + j];
                v_e[i + j] = v_c_e[i + j];
                w_e[i + j] = w_c_e[i + j];
            }
            col_count =0;
            i += 11;
        }
    }
    ///Readout loop:
    unsigned long count(0);
    unsigned long h_i(0);
    for(unsigned long i(6*(d_width+4)+6); i < 3*((d_height+2)*(d_width+4));)
    {
        if(count < d_width)
        {
            h_e[h_i] = u_e[i];
            ++h_i;
            ++count;
            ++i;
            ++i;
            ++i;
        }
        else
        {
            i += 12;
            count = 0;
        }
    }

}

    void CorrectionProcessing<boundaries::REFLECT, tags::CPU::SSE>
::value(DenseVector<double> & predictedu,
        DenseVector<double> & predictedv,
        DenseVector<double> & predictedw,
        DenseVector<double> & u,
        DenseVector<double> & v,
        DenseVector<double> & w,
        unsigned long d_width,
        unsigned long d_height,
        DenseMatrix<double> & height)
{
    double * u_e = u.elements();
    double * v_e = v.elements();
    double * w_e = w.elements();

    double * u_p_e = predictedu.elements();
    double * v_p_e = predictedv.elements();
    double * w_p_e = predictedw.elements();

    double * h_e = height.elements();

    unsigned long sse_steps((u.size() - (u.size() % 2)) / 2);
    unsigned long sse_limit(sse_steps * 2);
    double __attribute__((aligned(16))) v_c_e[v.size()];
    double __attribute__((aligned(16))) u_c_e[u.size()];
    double __attribute__((aligned(16))) w_c_e[w.size()];

    unsigned long c_count(0);
    for(unsigned long i(0); i < u.size(); ++i)
    {
        if(i < 6*(d_width+4)+6 || i >3 * ((d_height + 2) * (d_width + 4)))
        {

            v_c_e[i] = v_e[i];
            u_c_e[i] = u_e[i];
            w_c_e[i] = w_e[i];

            ++c_count;
        }
        else if(c_count == d_width)
        {
            for(unsigned j(0); j < 12; ++j)
            {
                u_c_e[i + j] = u_e[i + j];
                v_c_e[i + j] = v_e[i + j];
                w_c_e[i + j] = w_e[i + j];
            }
            c_count =0;
            i += 11;
        }

    }

    ///SSE loop:
    __m128d m1, m2, m3, m4, m5, m6, m7, m8;
    double __attribute__((aligned(16))) one_half(0.5);
    m1 = _mm_set_pd1(one_half);

    for(unsigned long i(0); i < sse_limit; i += 2)
    {
        m2 = _mm_load_pd(u_e + i);
        m3 = _mm_load_pd(v_e + i);
        m4 = _mm_load_pd(w_e + i);
        m5 = _mm_load_pd(u_p_e + i);
        m6 = _mm_load_pd(v_p_e + i);
        m7 = _mm_load_pd(w_p_e + i);

        m8 = _mm_add_pd(m2, m5);
        m2 = _mm_mul_pd(m8, m1);

        m8 = _mm_add_pd(m3, m6);
        m3 = _mm_mul_pd(m8, m1);

        m8 = _mm_add_pd(m4, m7);
        m4 = _mm_mul_pd(m8, m1);

        _mm_store_pd(u_e + i, m2);
        _mm_store_pd(v_e + i, m3);
        _mm_store_pd(w_e + i, m4);
    }

    ///Non SSE computation:
    for(unsigned long i(sse_limit); i < u.size(); ++i)
    {
        u_e[i] = double(0.5) * (u_e[i] + u_p_e[i]);
        v_e[i] = double(0.5) * (v_e[i] + v_p_e[i]);
        w_e[i] = double(0.5) * (w_e[i] + w_p_e[i]);
    }
    unsigned long col_count(0);
    for(unsigned long i(0); i < u.size(); ++i)
    {
        if(i < 6*(d_width+4)+6 || i >3 * ((d_height + 2) * (d_width + 4)))
        {
            u_e[i] = u_c_e[i];
            v_e[i] = v_c_e[i];
            w_e[i] = w_c_e[i];

            ++col_count;
        }
        else if(col_count == d_width)
        {
            for(unsigned j(0); j < 12; ++j)
            {
                u_e[i + j] = u_c_e[i + j];
                v_e[i + j] = v_c_e[i + j];
                w_e[i + j] = w_c_e[i + j];
            }
            col_count =0;
            i += 11;
        }
    }

    ///Readout loop:
    unsigned long count(0);
    unsigned long h_i(0);
    for(unsigned long i(6*(d_width+4)+6); i < 3*((d_height+2)*(d_width+4));)
    {
        if(count < d_width)
        {
            h_e[h_i] = u_e[i];
            ++h_i;
            ++count;
            ++i;
            ++i;
            ++i;
        }
        else
        {
            i += 12;
            count = 0;
        }
    }

}

