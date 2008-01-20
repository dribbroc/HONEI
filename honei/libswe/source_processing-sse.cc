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

#include <honei/libswe/source_processing.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

DenseVector<float> SourceProcessing<source_types::SIMPLE, tags::CPU::SSE>::value(DenseVector<float> & vector, DenseVector<float> & bottom_slopes_x, DenseVector<float> & bottom_slopes_y, float manning_n_squared)
{
    float * v_e = vector.elements();
    float * sl_x_e = bottom_slopes_x.elements();
    float * sl_y_e = bottom_slopes_y.elements();
    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 4)) / 4);
    unsigned long sse_limit(sse_tripels * 12);
    float __attribute__((aligned(16))) h_buffer[tripels];
    float __attribute__((aligned(16))) q1_buffer[tripels];
    float __attribute__((aligned(16))) q2_buffer[tripels];

    float __attribute__((aligned(16))) result_h[tripels];
    float __attribute__((aligned(16))) result_q1[tripels];
    float __attribute__((aligned(16))) result_q2[tripels];

    float __attribute__((aligned(16))) result_pow[tripels];

    float __attribute__((aligned(16))) g(9.81f);
    float __attribute__((aligned(16))) zero(0.f);
    float __attribute__((aligned(16))) exponent(float(-7.f)/float(3.f));
    float __attribute__((aligned(16))) minus_one(-1.f);
    float __attribute__((aligned(16))) manning(manning_n_squared);

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        //std::cout << v_e[i] << ", " << v_e[i + 1] << ", " << v_e[i+2] << endl;
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }
/*    ///Preprocess pow terms:
    buffer_i = 0;
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        float h = v_e[i];
        if (fabs(h) >= std::numeric_limits<float>::epsilon())
        {
            result_pow[buffer_i] = pow(h, exponent);
        }
        else
        {
            result_pow[buffer_i] = pow(std::numeric_limits<float>::epsilon(), exponent);
        }
        //std::cout << h << "^" << exponent << "=" <<result_pow[buffer_i] << endl;
        ++buffer_i;
    }
*/
    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) and uses manning = 0!!! TODO!!
    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    m1 = _mm_set_ps1(g);
    m2 = _mm_set_ps1(zero);
    m3 = _mm_set_ps1(exponent);
    m4 = _mm_set_ps1(minus_one);
    m5 = _mm_set_ps1(manning);

    for(unsigned long i(0); i < tripels - (tripels % 4); i += 4)
    {
        m6 = _mm_load_ps(h_buffer + i);
        m7 = _mm_load_ps(sl_x_e + i);
        m8 = _mm_load_ps(sl_y_e + i);

        m7 = _mm_mul_ps(m6, m7);
        m7 = _mm_mul_ps(m7, m4);
        m7 = _mm_mul_ps(m7, m1);

        m8 = _mm_mul_ps(m6, m8);
        m8 = _mm_mul_ps(m8, m4);
        m8 = _mm_mul_ps(m8, m1);

        _mm_store_ps(result_h + i, m2);
        _mm_store_ps(result_q1 + i, m7);
        _mm_store_ps(result_q2 + i, m8);
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
    buffer_i = 0;
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        float h = v_e[i];
        float q1 = v_e[i + 1];
        float q2 = v_e[i + 2];

        if (fabs(h) >= std::numeric_limits<float>::epsilon())
        {
            v_e[i] = float(0);
            //float temp = manning_n_squared * pow(h, float(float(-7)/float(3))) * sqrt(q1 * q1 + q2 * q2) * float(-1.);

            v_e[i + 1] = /*(temp * q1*/ - h * sl_x_e[buffer_i]/*)*/ * float(9.81);
            v_e[i + 2] = /*(temp * q2*/ - h * sl_y_e[buffer_i]/*)*/ * float(9.81);

            ++buffer_i;
        }
        else
        {
            v_e[i] = float(0);
            //float temp = manning_n_squared * pow(std::numeric_limits<float>::epsilon(), float(-7/3)) * sqrt(q1 * q1 + q2 * q2) * float(-1.);

            v_e[i + 1] = /*(temp * q1*/ - std::numeric_limits<float>::epsilon() * sl_x_e[buffer_i]/*)*/ * float(9.81);
            v_e[i + 2] = /*(temp * q2*/ - std::numeric_limits<float>::epsilon() * sl_y_e[buffer_i]/*)*/ * float(9.81);

            ++buffer_i;
        }
    }
    //std::cout << "n: " << vector.size() << endl;
    //std::cout << "Tripels: " << tripels << endl;
    //std::cout << "sse_limit: " << sse_limit << endl;
    return vector;
}
DenseVector<double> SourceProcessing<source_types::SIMPLE, tags::CPU::SSE>::value(DenseVector<double> & vector, DenseVector<double> & bottom_slopes_x, DenseVector<double> & bottom_slopes_y, double manning_n_squared)
{
    double * v_e = vector.elements();
    double * sl_x_e = bottom_slopes_x.elements();
    double * sl_y_e = bottom_slopes_y.elements();
    unsigned long tripels(vector.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 2)) / 2);
    unsigned long sse_limit(sse_tripels * 6);
    double __attribute__((aligned(16))) h_buffer[tripels];
    double __attribute__((aligned(16))) q1_buffer[tripels];
    double __attribute__((aligned(16))) q2_buffer[tripels];

    double __attribute__((aligned(16))) result_h[tripels];
    double __attribute__((aligned(16))) result_q1[tripels];
    double __attribute__((aligned(16))) result_q2[tripels];

    double __attribute__((aligned(16))) result_pow[tripels];

    double __attribute__((aligned(16))) g(9.81);
    double __attribute__((aligned(16))) zero(0.);
    double __attribute__((aligned(16))) exponent(double(-7.)/double(3.));
    double __attribute__((aligned(16))) minus_one(-1.);
    double __attribute__((aligned(16))) manning(manning_n_squared);

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        //std::cout << v_e[i] << ", " << v_e[i + 1] << ", " << v_e[i+2] << endl;
        h_buffer[buffer_i] = v_e[i];
        q1_buffer[buffer_i] = v_e[i + 1];
        q2_buffer[buffer_i] = v_e[i + 2];
        ++buffer_i;
    }
/*    ///Preprocess pow terms:
    buffer_i = 0;
    for(unsigned long i(0); i < vector.size(); i += 3)
    {
        float h = v_e[i];
        if (fabs(h) >= std::numeric_limits<float>::epsilon())
        {
            result_pow[buffer_i] = pow(h, exponent);
        }
        else
        {
            result_pow[buffer_i] = pow(std::numeric_limits<float>::epsilon(), exponent);
        }
        //std::cout << h << "^" << exponent << "=" <<result_pow[buffer_i] << endl;
        ++buffer_i;
    }
*/
    ///SSE loop: This routine does not yet check for any NAN cases (dry states!) and uses manning = 0!!! TODO!!
    __m128d m1, m2, m3, m4, m5, m6, m7, m8;
    m1 = _mm_set_pd1(g);
    m2 = _mm_set_pd1(zero);
    m3 = _mm_set_pd1(exponent);
    m4 = _mm_set_pd1(minus_one);
    m5 = _mm_set_pd1(manning);

    for(unsigned long i(0); i < tripels - (tripels % 2); i += 2)
    {
        m6 = _mm_load_pd(h_buffer + i);
        m7 = _mm_load_pd(sl_x_e + i);
        m8 = _mm_load_pd(sl_y_e + i);

        m7 = _mm_mul_pd(m6, m7);
        m7 = _mm_mul_pd(m7, m4);
        m7 = _mm_mul_pd(m7, m1);

        m8 = _mm_mul_pd(m6, m8);
        m8 = _mm_mul_pd(m8, m4);
        m8 = _mm_mul_pd(m8, m1);

        _mm_store_pd(result_h + i, m2);
        _mm_store_pd(result_q1 + i, m7);
        _mm_store_pd(result_q2 + i, m8);
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
    buffer_i = 0;
    for(unsigned long i(sse_limit); i < vector.size(); i += 3)
    {
        double h = v_e[i];
        double q1 = v_e[i + 1];
        double q2 = v_e[i + 2];

        if (fabs(h) >= std::numeric_limits<float>::epsilon())
        {
            v_e[i] = double(0);
            //float temp = manning_n_squared * pow(h, float(float(-7)/float(3))) * sqrt(q1 * q1 + q2 * q2) * float(-1.);

            v_e[i + 1] = /*(temp * q1*/ - h * sl_x_e[buffer_i]/*)*/ * double(9.81);
            v_e[i + 2] = /*(temp * q2*/ - h * sl_y_e[buffer_i]/*)*/ * double(9.81);

            ++buffer_i;
        }
        else
        {
            v_e[i] = double(0);
            //float temp = manning_n_squared * pow(std::numeric_limits<float>::epsilon(), float(-7/3)) * sqrt(q1 * q1 + q2 * q2) * float(-1.);

            v_e[i + 1] = /*(temp * q1*/ - std::numeric_limits<double>::epsilon() * sl_x_e[buffer_i]/*)*/ * double(9.81);
            v_e[i + 2] = /*(temp * q2*/ - std::numeric_limits<double>::epsilon() * sl_y_e[buffer_i]/*)*/ * double(9.81);

            ++buffer_i;
        }
    }
    //std::cout << "n: " << vector.size() << endl;
    //std::cout << "Tripels: " << tripels << endl;
    //std::cout << "sse_limit: " << sse_limit << endl;
    return vector;
}

