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

#include <honei/util/attributes.hh>
#include <honei/libswe/flow_processing.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

DenseVector<float> FlowProcessing<directions::X, tags::CPU::SSE>::value(DenseVector<float> & vector)
{
    float * v_e = vector.elements();

    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 4)) / 4);
    unsigned long sse_limit(sse_tripels * 12);
    float HONEI_ALIGNED(16) h_buffer[tripels];
    float HONEI_ALIGNED(16) q1_buffer[tripels];
    float HONEI_ALIGNED(16) q2_buffer[tripels];

    float HONEI_ALIGNED(16) result_h[tripels];
    float HONEI_ALIGNED(16) result_q1[tripels];
    float HONEI_ALIGNED(16) result_q2[tripels];

    ///Alignment loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }

    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) TODO!!
    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float HONEI_ALIGNED(16) g(9.81f);
    float HONEI_ALIGNED(16) two(2.0f);
    m4 = _mm_set_ps1(g);
    m5 = _mm_set_ps1(two);
    for(unsigned long i(0); i < tripels - (tripels % 4); i += 4)
    {
        m1 = _mm_load_ps(h_buffer + i);
        m2 = _mm_load_ps(q1_buffer + i);
        m3 = _mm_load_ps(q2_buffer + i);

        // a = q1 * q1 / h
        m6 = _mm_mul_ps(m2, m2);
        m6 = _mm_div_ps(m6, m1);

        // b = g * h^2 / 2
        m7 = _mm_mul_ps(m1, m1);
        m7 = _mm_mul_ps(m7, m4);
        m7 = _mm_div_ps(m7, m5);

        //a + b
        m6 = _mm_add_ps(m7, m6);

        //q1 * q2 / h
        m7 = _mm_mul_ps(m2, m3);
        m7 = _mm_div_ps(m7, m1);

        ///TODO: stream?
        _mm_store_ps(result_h + i, m2);
        _mm_store_ps(result_q1 + i, m6);
        _mm_store_ps(result_q2 + i, m7);
    }

    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        v_e[i] = result_h[buffer_i];
        v_e[i + 1] = result_q1[buffer_i];
        v_e[i + 2] = result_q2[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        if (fabs(v_e[i]) >= std::numeric_limits<float>::epsilon())
        {

            float h = v_e[i];
            float q1 = v_e[i + 1];
            float q2 = v_e[i + 2];

            v_e[i] = q1;
            v_e[i + 1] = (q1 * q1 / h) + (float(9.81) * h * h / float(2));
            v_e[i + 2] = q1 * q2 / h;
        }
        else
        {
            float h = v_e[i];
            float q1 = v_e[i + 1];
            float q2 = v_e[i + 2];

            v_e[i] = q1;
            v_e[i + 1] = q1 * q1 / std::numeric_limits<float>::epsilon();
            v_e[i + 2] = q1 * q2/ std::numeric_limits<float>::epsilon();
        }
    }

    return vector;
}

DenseVector<float> FlowProcessing<directions::Y, tags::CPU::SSE>::value(DenseVector<float> & vector)
{
    float * v_e = vector.elements();

    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 4)) / 4);
    unsigned long sse_limit(sse_tripels * 12);
    float HONEI_ALIGNED(16) h_buffer[tripels];
    float HONEI_ALIGNED(16) q1_buffer[tripels];
    float HONEI_ALIGNED(16) q2_buffer[tripels];

    float HONEI_ALIGNED(16) result_h[tripels];
    float HONEI_ALIGNED(16) result_q1[tripels];
    float HONEI_ALIGNED(16) result_q2[tripels];

    ///Alignment loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }

    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) TODO!!
    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float HONEI_ALIGNED(16) g(9.81f);
    float HONEI_ALIGNED(16) two(2.0f);
    m4 = _mm_set_ps1(g);
    m5 = _mm_set_ps1(two);
    for(unsigned long i(0); i < tripels - (tripels % 4); i += 4)
    {
        m1 = _mm_load_ps(h_buffer + i);
        m2 = _mm_load_ps(q1_buffer + i);
        m3 = _mm_load_ps(q2_buffer + i);

        // a = q2 * q2 / h
        m6 = _mm_mul_ps(m3, m3);
        m6 = _mm_div_ps(m6, m1);

        // b = g * h^2 / 2
        m7 = _mm_mul_ps(m1, m1);
        m7 = _mm_mul_ps(m7, m4);
        m7 = _mm_div_ps(m7, m5);

        //a + b
        m6 = _mm_add_ps(m7, m6);

        //q1 * q2 / h
        m7 = _mm_mul_ps(m2, m3);
        m7 = _mm_div_ps(m7, m1);

        ///TODO: stream?
        _mm_store_ps(result_h + i, m3);
        _mm_store_ps(result_q1 + i, m7);
        _mm_store_ps(result_q2 + i, m6);
    }

    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        v_e[i] = result_h[buffer_i];
        v_e[i + 1] = result_q1[buffer_i];
        v_e[i + 2] = result_q2[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        if (fabs(v_e[i]) >= std::numeric_limits<float>::epsilon())
        {
            float h = v_e[i];
            float q1 = v_e[i + 1];
            float q2 = v_e[i + 2];

            v_e[i] = q2;
            v_e[i + 1] = q1 * q2 / h ;
            v_e[i + 2] = (q2 * q2 / h) + (float(9.81) * h * h / float(2));
        }
        else
        {
            float h = v_e[i];
            float q1 = v_e[i + 1];
            float q2 = v_e[i + 2];

            v_e[i] = q2;
            v_e[i + 1] = q1 * q2 / h ;
            v_e[i + 2] = (q2 * q2 / std::numeric_limits<float>::epsilon());
        }
    }

    return vector;
}

DenseVector<double> FlowProcessing<directions::X, tags::CPU::SSE>::value(DenseVector<double> & vector)
{
    double * v_e = vector.elements();

    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 2)) / 2);
    unsigned long sse_limit(sse_tripels * 6);
    double HONEI_ALIGNED(16) h_buffer[tripels];
    double HONEI_ALIGNED(16) q1_buffer[tripels];
    double HONEI_ALIGNED(16) q2_buffer[tripels];

    double HONEI_ALIGNED(16) result_h[tripels];
    double HONEI_ALIGNED(16) result_q1[tripels];
    double HONEI_ALIGNED(16) result_q2[tripels];

    ///Alignment loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }

    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) TODO!!
    __m128d m1, m2, m3, m4, m5, m6, m7, m8;
    double HONEI_ALIGNED(16) g(9.81);
    double HONEI_ALIGNED(16) two(2.0);
    m4 = _mm_set_pd1(g);
    m5 = _mm_set_pd1(two);
    for(unsigned long i(0); i < tripels - (tripels % 4); i += 2)
    {
        m1 = _mm_load_pd(h_buffer + i);
        m2 = _mm_load_pd(q1_buffer + i);
        m3 = _mm_load_pd(q2_buffer + i);

        // a = q1 * q1 / h
        m6 = _mm_mul_pd(m2, m2);
        m6 = _mm_div_pd(m6, m1);

        // b = g * h^2 / 2
        m7 = _mm_mul_pd(m1, m1);
        m7 = _mm_mul_pd(m7, m4);
        m7 = _mm_div_pd(m7, m5);

        //a + b
        m6 = _mm_add_pd(m7, m6);

        //q1 * q2 / h
        m7 = _mm_mul_pd(m2, m3);
        m7 = _mm_div_pd(m7, m1);

        ///TODO: stream?
        _mm_store_pd(result_h + i, m2);
        _mm_store_pd(result_q1 + i, m6);
        _mm_store_pd(result_q2 + i, m7);
    }

    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        v_e[i] = result_h[buffer_i];
        v_e[i + 1] = result_q1[buffer_i];
        v_e[i + 2] = result_q2[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        if (fabs(v_e[i]) >= std::numeric_limits<double>::epsilon())
        {

            double h = v_e[i];
            double q1 = v_e[i + 1];
            double q2 = v_e[i + 2];

            v_e[i] = q1;
            v_e[i + 1] = (q1 * q1 / h) + (double(9.81) * h * h / double(2));
            v_e[i + 2] = q1 * q2 / h;
        }
        else
        {
            double h = v_e[i];
            double q1 = v_e[i + 1];
            double q2 = v_e[i + 2];

            v_e[i] = q1;
            v_e[i + 1] = q1 * q1 / std::numeric_limits<double>::epsilon();
            v_e[i + 2] = q1 * q2/ std::numeric_limits<double>::epsilon();
        }
    }

    return vector;
}

DenseVector<double> FlowProcessing<directions::Y, tags::CPU::SSE>::value(DenseVector<double> & vector)
{
    double * v_e = vector.elements();

    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 2)) / 2);
    unsigned long sse_limit(sse_tripels * 6);
    double HONEI_ALIGNED(16) h_buffer[tripels];
    double HONEI_ALIGNED(16) q1_buffer[tripels];
    double HONEI_ALIGNED(16) q2_buffer[tripels];

    double HONEI_ALIGNED(16) result_h[tripels];
    double HONEI_ALIGNED(16) result_q1[tripels];
    double HONEI_ALIGNED(16) result_q2[tripels];

    ///Alignment loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }

    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) TODO!!
    __m128d m1, m2, m3, m4, m5, m6, m7, m8;
    double HONEI_ALIGNED(16) g(9.81);
    double HONEI_ALIGNED(16) two(2.0);
    m4 = _mm_set_pd1(g);
    m5 = _mm_set_pd1(two);
    for(unsigned long i(0); i < tripels - (tripels % 4); i += 2)
    {
        m1 = _mm_load_pd(h_buffer + i);
        m2 = _mm_load_pd(q1_buffer + i);
        m3 = _mm_load_pd(q2_buffer + i);

        // a = q2 * q2 / h
        m6 = _mm_mul_pd(m3, m3);
        m6 = _mm_div_pd(m6, m1);

        // b = g * h^2 / 2
        m7 = _mm_mul_pd(m1, m1);
        m7 = _mm_mul_pd(m7, m4);
        m7 = _mm_div_pd(m7, m5);

        //a + b
        m6 = _mm_add_pd(m7, m6);

        //q1 * q2 / h
        m7 = _mm_mul_pd(m2, m3);
        m7 = _mm_div_pd(m7, m1);

        ///TODO: stream?
        _mm_store_pd(result_h + i, m3);
        _mm_store_pd(result_q1 + i, m7);
        _mm_store_pd(result_q2 + i, m6);
    }

    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        v_e[i] = result_h[buffer_i];
        v_e[i + 1] = result_q1[buffer_i];
        v_e[i + 2] = result_q2[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        if (fabs(v_e[i]) >= std::numeric_limits<double>::epsilon())
        {
            double h = v_e[i];
            double q1 = v_e[i + 1];
            double q2 = v_e[i + 2];

            v_e[i] = q2;
            v_e[i + 1] = q1 * q2 / h ;
            v_e[i + 2] = (q2 * q2 / h) + (double(9.81) * h * h / double(2));
        }
        else
        {
            double h = v_e[i];
            double q1 = v_e[i + 1];
            double q2 = v_e[i + 2];

            v_e[i] = q2;
            v_e[i + 1] = q1 * q2 / h ;
            v_e[i + 2] = (q2 * q2 / std::numeric_limits<double>::epsilon());
        }
    }

    return vector;
}

