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

#include <honei/attributes.hh>
#include <honei/libswe/assembly_processing.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

BandedMatrix<float> & AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M6>::value(BandedMatrix<float> & m1, BandedMatrix<float> & result, DenseVector<float> & c, unsigned long d_width, unsigned long d_height)
{
    DenseVector<float> m6diag(m1.size());
    DenseVector<float> m6bandplus3(m1.size());
    DenseVector<float> m6bandplus6(m1.size());
    DenseVector<float> m6bandminus3(m1.size());

    float * m6diag_e = m6diag.elements();
    float * m6bandplus3_e = m6bandplus3.elements();
    float * m6bandplus6_e = m6bandplus6.elements();
    float * m6bandminus3_e = m6bandminus3.elements();

    float * m1diag_e = (m1.band(0)).elements();
    float * m1bandplus3_e = (m1.band(3)).elements();
    float * m1bandplus6_e = (m1.band(6)).elements();
    float * m1bandminus3_e = (m1.band(-3)).elements();

    unsigned long tripels(m1.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 4)) / 4);
    unsigned long sse_limit(sse_tripels * 12);

    float HONEI_ATTRIBUTE(aligned(16)) c1_buffer[4];
    float HONEI_ATTRIBUTE(aligned(16)) c2_buffer[4];
    float HONEI_ATTRIBUTE(aligned(16)) c3_buffer[4];

    for(int i(0); i < 4; ++i)
    {
        c1_buffer[i] = c[0];
        c2_buffer[i] = c[1];
        c3_buffer[i] = c[2];
    }

    float HONEI_ATTRIBUTE(aligned(16)) d_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) d_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) d_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) b1_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b1_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b1_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) b2_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b2_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b2_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) bm1_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) bm1_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) bm1_q2_buffer[tripels];

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < m1.size(); i += 3)
    {
        d_h_buffer[buffer_i] = m1diag_e[i];
        d_q1_buffer[buffer_i] = m1diag_e[i + 1];
        d_q2_buffer[buffer_i] = m1diag_e[i + 2];
        b1_h_buffer[buffer_i] = m1bandplus3_e[i];
        b1_q1_buffer[buffer_i] = m1bandplus3_e[i + 1];
        b1_q2_buffer[buffer_i] = m1bandplus3_e[i + 2];
        b2_h_buffer[buffer_i] = m1bandplus6_e[i];
        b2_q1_buffer[buffer_i] = m1bandplus6_e[i + 1];
        b2_q2_buffer[buffer_i] = m1bandplus6_e[i + 2];
        bm1_h_buffer[buffer_i] = m1bandminus3_e[i];
        bm1_q1_buffer[buffer_i] = m1bandminus3_e[i + 1];
        bm1_q2_buffer[buffer_i] = m1bandminus3_e[i + 2];
        ++buffer_i;
    }

    __m128 m0, m2, m3, m4, m5, m6, m7, m8;

    m0 = _mm_load_ps(c1_buffer);
    m2 = _mm_load_ps(c2_buffer);
    m3 = _mm_load_ps(c3_buffer);

    m0 = _mm_mul_ps(m0, m0);
    m2 = _mm_mul_ps(m2, m2);
    m3 = _mm_mul_ps(m3, m3);

    _mm_store_ps(c1_buffer, m0);
    _mm_store_ps(c2_buffer, m2);
    _mm_store_ps(c3_buffer, m3);


    for(unsigned long i(0); i < tripels - (tripels % 4); i += 4)
    {
        m4 = _mm_load_ps(d_h_buffer + i);
        m5 = _mm_load_ps(b1_h_buffer + i);
        m6 = _mm_load_ps(b2_h_buffer + i);
        m7 = _mm_load_ps(bm1_h_buffer + i);

        m4 = _mm_mul_ps(m4, m0);
        m5 = _mm_mul_ps(m5, m0);
        m6 = _mm_mul_ps(m6, m0);
        m7 = _mm_mul_ps(m7, m0);

        _mm_store_ps(d_h_buffer + i, m4);
        _mm_store_ps(b1_h_buffer + i, m5);
        _mm_store_ps(b2_h_buffer + i, m6);
        _mm_store_ps(bm1_h_buffer + i, m7);

        m4 = _mm_load_ps(d_q1_buffer + i);
        m5 = _mm_load_ps(b1_q1_buffer + i);
        m6 = _mm_load_ps(b2_q1_buffer + i);
        m7 = _mm_load_ps(bm1_q1_buffer + i);

        m4 = _mm_mul_ps(m4, m2);
        m5 = _mm_mul_ps(m5, m2);
        m6 = _mm_mul_ps(m6, m2);
        m7 = _mm_mul_ps(m7, m2);

        _mm_store_ps(d_q1_buffer + i, m4);
        _mm_store_ps(b1_q1_buffer + i, m5);
        _mm_store_ps(b2_q1_buffer + i, m6);
        _mm_store_ps(bm1_q1_buffer + i, m7);

        m4 = _mm_load_ps(d_q2_buffer + i);
        m5 = _mm_load_ps(b1_q2_buffer + i);
        m6 = _mm_load_ps(b2_q2_buffer + i);
        m7 = _mm_load_ps(bm1_q2_buffer + i);

        m4 = _mm_mul_ps(m4, m3);
        m5 = _mm_mul_ps(m5, m3);
        m6 = _mm_mul_ps(m6, m3);
        m7 = _mm_mul_ps(m7, m3);

        _mm_store_ps(d_q2_buffer + i, m4);
        _mm_store_ps(b1_q2_buffer + i, m5);
        _mm_store_ps(b2_q2_buffer + i, m6);
        _mm_store_ps(bm1_q2_buffer + i, m7);
    }
    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        m6diag_e[i] = d_h_buffer[buffer_i];
        m6diag_e[i + 1] = d_q1_buffer[buffer_i];
        m6diag_e[i + 2] = d_q2_buffer[buffer_i];
        m6bandplus3_e[i] = b1_h_buffer[buffer_i];
        m6bandplus3_e[i + 1] = b1_q1_buffer[buffer_i];
        m6bandplus3_e[i + 2] = b1_q2_buffer[buffer_i];
        m6bandplus6_e[i] = b2_h_buffer[buffer_i];
        m6bandplus6_e[i + 1] = b2_q1_buffer[buffer_i];
        m6bandplus6_e[i + 2] = b2_q2_buffer[buffer_i];
        m6bandminus3_e[i] = bm1_h_buffer[buffer_i];
        m6bandminus3_e[i + 1] = bm1_q1_buffer[buffer_i];
        m6bandminus3_e[i + 2] = bm1_q2_buffer[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < m1.size(); i += 3)
    {
        m6diag_e[i] = c1_buffer[0] * m1diag_e[i];
        m6diag_e[i + 1] = c2_buffer[0] * m1diag_e[i + 1];
        m6diag_e[i + 2] = c3_buffer[0] * m1diag_e[i + 2];
        m6bandplus3_e[i] = c1_buffer[0] * m1bandplus3_e[i];
        m6bandplus3_e[i + 1] = c2_buffer[0] * m1bandplus3_e[i + 1];
        m6bandplus3_e[i + 2] = c3_buffer[0] * m1bandplus3_e[i + 2];
        m6bandplus6_e[i] = c1_buffer[0] * m1bandplus6_e[i];
        m6bandplus6_e[i + 1] = c2_buffer[0] * m1bandplus6_e[i + 1];
        m6bandplus6_e[i + 2] = c3_buffer[0] * m1bandplus6_e[i + 2];
        m6bandminus3_e[i] = c1_buffer[0] * m1bandminus3_e[i];
        m6bandminus3_e[i + 1] = c2_buffer[0] * m1bandminus3_e[i + 1];
        m6bandminus3_e[i + 2] = c3_buffer[0] * m1bandminus3_e[i + 2];
    }

    result.insert_band(0, m6diag);
    result.insert_band(3, m6bandplus3);
    result.insert_band(6, m6bandplus6);
    result.insert_band(-3, m6bandminus3);

    return result;
}

BandedMatrix<double> & AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M6>::value(BandedMatrix<double> & m1, BandedMatrix<double> & result, DenseVector<double> & c, unsigned long d_width, unsigned long d_height)
{
    DenseVector<double> m6diag(m1.size());
    DenseVector<double> m6bandplus3(m1.size());
    DenseVector<double> m6bandplus6(m1.size());
    DenseVector<double> m6bandminus3(m1.size());

    double * m6diag_e = m6diag.elements();
    double * m6bandplus3_e = m6bandplus3.elements();
    double * m6bandplus6_e = m6bandplus6.elements();
    double * m6bandminus3_e = m6bandminus3.elements();

    double * m1diag_e = (m1.band(0)).elements();
    double * m1bandplus3_e = (m1.band(3)).elements();
    double * m1bandplus6_e = (m1.band(6)).elements();
    double * m1bandminus3_e = (m1.band(-3)).elements();

    unsigned long tripels(m1.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 2)) / 2);
    unsigned long sse_limit(sse_tripels * 6);

    double HONEI_ATTRIBUTE(aligned(16)) c1_buffer[4];
    double HONEI_ATTRIBUTE(aligned(16)) c2_buffer[4];
    double HONEI_ATTRIBUTE(aligned(16)) c3_buffer[4];

    for(int i(0); i < 2; ++i)
    {
        c1_buffer[i] = c[0];
        c2_buffer[i] = c[1];
        c3_buffer[i] = c[2];
    }

    double HONEI_ATTRIBUTE(aligned(16)) d_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) d_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) d_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) b1_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b1_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b1_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) b2_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b2_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b2_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) bm1_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) bm1_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) bm1_q2_buffer[tripels];

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < m1.size(); i += 3)
    {
        d_h_buffer[buffer_i] = m1diag_e[i];
        d_q1_buffer[buffer_i] = m1diag_e[i + 1];
        d_q2_buffer[buffer_i] = m1diag_e[i + 2];
        b1_h_buffer[buffer_i] = m1bandplus3_e[i];
        b1_q1_buffer[buffer_i] = m1bandplus3_e[i + 1];
        b1_q2_buffer[buffer_i] = m1bandplus3_e[i + 2];
        b2_h_buffer[buffer_i] = m1bandplus6_e[i];
        b2_q1_buffer[buffer_i] = m1bandplus6_e[i + 1];
        b2_q2_buffer[buffer_i] = m1bandplus6_e[i + 2];
        bm1_h_buffer[buffer_i] = m1bandminus3_e[i];
        bm1_q1_buffer[buffer_i] = m1bandminus3_e[i + 1];
        bm1_q2_buffer[buffer_i] = m1bandminus3_e[i + 2];
        ++buffer_i;
    }

    __m128d m0, m2, m3, m4, m5, m6, m7, m8;

    m0 = _mm_load_pd(c1_buffer);
    m2 = _mm_load_pd(c2_buffer);
    m3 = _mm_load_pd(c3_buffer);

    m0 = _mm_mul_pd(m0, m0);
    m2 = _mm_mul_pd(m2, m2);
    m3 = _mm_mul_pd(m3, m3);

    _mm_store_pd(c1_buffer, m0);
    _mm_store_pd(c2_buffer, m2);
    _mm_store_pd(c3_buffer, m3);


    for(unsigned long i(0); i < tripels - (tripels % 2); i += 2)
    {
        m4 = _mm_load_pd(d_h_buffer + i);
        m5 = _mm_load_pd(b1_h_buffer + i);
        m6 = _mm_load_pd(b2_h_buffer + i);
        m7 = _mm_load_pd(bm1_h_buffer + i);

        m4 = _mm_mul_pd(m4, m0);
        m5 = _mm_mul_pd(m5, m0);
        m6 = _mm_mul_pd(m6, m0);
        m7 = _mm_mul_pd(m7, m0);

        _mm_store_pd(d_h_buffer + i, m4);
        _mm_store_pd(b1_h_buffer + i, m5);
        _mm_store_pd(b2_h_buffer + i, m6);
        _mm_store_pd(bm1_h_buffer + i, m7);

        m4 = _mm_load_pd(d_q1_buffer + i);
        m5 = _mm_load_pd(b1_q1_buffer + i);
        m6 = _mm_load_pd(b2_q1_buffer + i);
        m7 = _mm_load_pd(bm1_q1_buffer + i);

        m4 = _mm_mul_pd(m4, m2);
        m5 = _mm_mul_pd(m5, m2);
        m6 = _mm_mul_pd(m6, m2);
        m7 = _mm_mul_pd(m7, m2);

        _mm_store_pd(d_q1_buffer + i, m4);
        _mm_store_pd(b1_q1_buffer + i, m5);
        _mm_store_pd(b2_q1_buffer + i, m6);
        _mm_store_pd(bm1_q1_buffer + i, m7);

        m4 = _mm_load_pd(d_q2_buffer + i);
        m5 = _mm_load_pd(b1_q2_buffer + i);
        m6 = _mm_load_pd(b2_q2_buffer + i);
        m7 = _mm_load_pd(bm1_q2_buffer + i);

        m4 = _mm_mul_pd(m4, m3);
        m5 = _mm_mul_pd(m5, m3);
        m6 = _mm_mul_pd(m6, m3);
        m7 = _mm_mul_pd(m7, m3);

        _mm_store_pd(d_q2_buffer + i, m4);
        _mm_store_pd(b1_q2_buffer + i, m5);
        _mm_store_pd(b2_q2_buffer + i, m6);
        _mm_store_pd(bm1_q2_buffer + i, m7);
    }
    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        m6diag_e[i] = d_h_buffer[buffer_i];
        m6diag_e[i + 1] = d_q1_buffer[buffer_i];
        m6diag_e[i + 2] = d_q2_buffer[buffer_i];
        m6bandplus3_e[i] = b1_h_buffer[buffer_i];
        m6bandplus3_e[i + 1] = b1_q1_buffer[buffer_i];
        m6bandplus3_e[i + 2] = b1_q2_buffer[buffer_i];
        m6bandplus6_e[i] = b2_h_buffer[buffer_i];
        m6bandplus6_e[i + 1] = b2_q1_buffer[buffer_i];
        m6bandplus6_e[i + 2] = b2_q2_buffer[buffer_i];
        m6bandminus3_e[i] = bm1_h_buffer[buffer_i];
        m6bandminus3_e[i + 1] = bm1_q1_buffer[buffer_i];
        m6bandminus3_e[i + 2] = bm1_q2_buffer[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < m1.size(); i += 3)
    {
        m6diag_e[i] = c1_buffer[0] * m1diag_e[i];
        m6diag_e[i + 1] = c2_buffer[0] * m1diag_e[i + 1];
        m6diag_e[i + 2] = c3_buffer[0] * m1diag_e[i + 2];
        m6bandplus3_e[i] = c1_buffer[0] * m1bandplus3_e[i];
        m6bandplus3_e[i + 1] = c2_buffer[0] * m1bandplus3_e[i + 1];
        m6bandplus3_e[i + 2] = c3_buffer[0] * m1bandplus3_e[i + 2];
        m6bandplus6_e[i] = c1_buffer[0] * m1bandplus6_e[i];
        m6bandplus6_e[i + 1] = c2_buffer[0] * m1bandplus6_e[i + 1];
        m6bandplus6_e[i + 2] = c3_buffer[0] * m1bandplus6_e[i + 2];
        m6bandminus3_e[i] = c1_buffer[0] * m1bandminus3_e[i];
        m6bandminus3_e[i + 1] = c2_buffer[0] * m1bandminus3_e[i + 1];
        m6bandminus3_e[i + 2] = c3_buffer[0] * m1bandminus3_e[i + 2];
    }

    result.insert_band(0, m6diag);
    result.insert_band(3, m6bandplus3);
    result.insert_band(6, m6bandplus6);
    result.insert_band(-3, m6bandminus3);

    return result;
}

//------------M8------------

BandedMatrix<float> & AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M8>::value(BandedMatrix<float> & m2, BandedMatrix<float> & result, DenseVector<float> & d, unsigned long d_width, unsigned long d_height)
{
    DenseVector<float> m8diag(m2.size());
    DenseVector<float> m8bandplus3(m2.size());
    DenseVector<float> m8bandplus6(m2.size());
    DenseVector<float> m8bandminus3(m2.size());

    float * m8diag_e = m8diag.elements();
    float * m8bandplus3_e = m8bandplus3.elements();
    float * m8bandplus6_e = m8bandplus6.elements();
    float * m8bandminus3_e = m8bandminus3.elements();

    float * m2diag_e = (m2.band(0)).elements();
    float * m2bandplus3_e = (m2.band(3)).elements();
    float * m2bandplus6_e = (m2.band(6)).elements();
    float * m2bandminus3_e = (m2.band(-3)).elements();

    unsigned long tripels(m2.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 4)) / 4);
    unsigned long sse_limit(sse_tripels * 12);

    float HONEI_ATTRIBUTE(aligned(16)) d1_buffer[4];
    float HONEI_ATTRIBUTE(aligned(16)) d2_buffer[4];
    float HONEI_ATTRIBUTE(aligned(16)) d3_buffer[4];

    for(int i(0); i < 4; ++i)
    {
        d1_buffer[i] = d[0];
        d2_buffer[i] = d[1];
        d3_buffer[i] = d[2];
    }

    float HONEI_ATTRIBUTE(aligned(16)) d_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) d_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) d_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) b1_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b1_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b1_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) b2_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b2_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) b2_q2_buffer[tripels];

    float HONEI_ATTRIBUTE(aligned(16)) bm1_h_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) bm1_q1_buffer[tripels];
    float HONEI_ATTRIBUTE(aligned(16)) bm1_q2_buffer[tripels];

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < m2.size(); i += 3)
    {
        d_h_buffer[buffer_i] = m2diag_e[i];
        d_q1_buffer[buffer_i] = m2diag_e[i + 1];
        d_q2_buffer[buffer_i] = m2diag_e[i + 2];
        b1_h_buffer[buffer_i] = m2bandplus3_e[i];
        b1_q1_buffer[buffer_i] = m2bandplus3_e[i + 1];
        b1_q2_buffer[buffer_i] = m2bandplus3_e[i + 2];
        b2_h_buffer[buffer_i] = m2bandplus6_e[i];
        b2_q1_buffer[buffer_i] = m2bandplus6_e[i + 1];
        b2_q2_buffer[buffer_i] = m2bandplus6_e[i + 2];
        bm1_h_buffer[buffer_i] = m2bandminus3_e[i];
        bm1_q1_buffer[buffer_i] = m2bandminus3_e[i + 1];
        bm1_q2_buffer[buffer_i] = m2bandminus3_e[i + 2];
        ++buffer_i;
    }

    __m128 m0, m1, m3, m4, m5, m6, m7, m8;

    m0 = _mm_load_ps(d1_buffer);
    m1 = _mm_load_ps(d2_buffer);
    m3 = _mm_load_ps(d3_buffer);

    m0 = _mm_mul_ps(m0, m0);
    m1 = _mm_mul_ps(m1, m1);
    m3 = _mm_mul_ps(m3, m3);

    _mm_store_ps(d1_buffer, m0);
    _mm_store_ps(d2_buffer, m1);
    _mm_store_ps(d3_buffer, m3);


    for(unsigned long i(0); i < tripels - (tripels % 4); i += 4)
    {
        m4 = _mm_load_ps(d_h_buffer + i);
        m5 = _mm_load_ps(b1_h_buffer + i);
        m6 = _mm_load_ps(b2_h_buffer + i);
        m7 = _mm_load_ps(bm1_h_buffer + i);

        m4 = _mm_mul_ps(m4, m0);
        m5 = _mm_mul_ps(m5, m0);
        m6 = _mm_mul_ps(m6, m0);
        m7 = _mm_mul_ps(m7, m0);

        _mm_store_ps(d_h_buffer + i, m4);
        _mm_store_ps(b1_h_buffer + i, m5);
        _mm_store_ps(b2_h_buffer + i, m6);
        _mm_store_ps(bm1_h_buffer + i, m7);

        m4 = _mm_load_ps(d_q1_buffer + i);
        m5 = _mm_load_ps(b1_q1_buffer + i);
        m6 = _mm_load_ps(b2_q1_buffer + i);
        m7 = _mm_load_ps(bm1_q1_buffer + i);

        m4 = _mm_mul_ps(m4, m1);
        m5 = _mm_mul_ps(m5, m1);
        m6 = _mm_mul_ps(m6, m1);
        m7 = _mm_mul_ps(m7, m1);

        _mm_store_ps(d_q1_buffer + i, m4);
        _mm_store_ps(b1_q1_buffer + i, m5);
        _mm_store_ps(b2_q1_buffer + i, m6);
        _mm_store_ps(bm1_q1_buffer + i, m7);

        m4 = _mm_load_ps(d_q2_buffer + i);
        m5 = _mm_load_ps(b1_q2_buffer + i);
        m6 = _mm_load_ps(b2_q2_buffer + i);
        m7 = _mm_load_ps(bm1_q2_buffer + i);

        m4 = _mm_mul_ps(m4, m3);
        m5 = _mm_mul_ps(m5, m3);
        m6 = _mm_mul_ps(m6, m3);
        m7 = _mm_mul_ps(m7, m3);

        _mm_store_ps(d_q2_buffer + i, m4);
        _mm_store_ps(b1_q2_buffer + i, m5);
        _mm_store_ps(b2_q2_buffer + i, m6);
        _mm_store_ps(bm1_q2_buffer + i, m7);
    }
    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        m8diag_e[i] = d_h_buffer[buffer_i];
        m8diag_e[i + 1] = d_q1_buffer[buffer_i];
        m8diag_e[i + 2] = d_q2_buffer[buffer_i];
        m8bandplus3_e[i] = b1_h_buffer[buffer_i];
        m8bandplus3_e[i + 1] = b1_q1_buffer[buffer_i];
        m8bandplus3_e[i + 2] = b1_q2_buffer[buffer_i];
        m8bandplus6_e[i] = b2_h_buffer[buffer_i];
        m8bandplus6_e[i + 1] = b2_q1_buffer[buffer_i];
        m8bandplus6_e[i + 2] = b2_q2_buffer[buffer_i];
        m8bandminus3_e[i] = bm1_h_buffer[buffer_i];
        m8bandminus3_e[i + 1] = bm1_q1_buffer[buffer_i];
        m8bandminus3_e[i + 2] = bm1_q2_buffer[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < m2.size(); i += 3)
    {
        m8diag_e[i] = d1_buffer[0] * m2diag_e[i];
        m8diag_e[i + 1] = d2_buffer[0] * m2diag_e[i + 1];
        m8diag_e[i + 2] = d3_buffer[0] * m2diag_e[i + 2];
        m8bandplus3_e[i] = d1_buffer[0] * m2bandplus3_e[i];
        m8bandplus3_e[i + 1] = d2_buffer[0] * m2bandplus3_e[i + 1];
        m8bandplus3_e[i + 2] = d3_buffer[0] * m2bandplus3_e[i + 2];
        m8bandplus6_e[i] = d1_buffer[0] * m2bandplus6_e[i];
        m8bandplus6_e[i + 1] = d2_buffer[0] * m2bandplus6_e[i + 1];
        m8bandplus6_e[i + 2] = d3_buffer[0] * m2bandplus6_e[i + 2];
        m8bandminus3_e[i] = d1_buffer[0] * m2bandminus3_e[i];
        m8bandminus3_e[i + 1] = d2_buffer[0] * m2bandminus3_e[i + 1];
        m8bandminus3_e[i + 2] = d3_buffer[0] * m2bandminus3_e[i + 2];
    }

    result.insert_band(0, m8diag);
    result.insert_band(3, m8bandplus3);
    result.insert_band(6, m8bandplus6);
    result.insert_band(-3, m8bandminus3);

    return result;
}

BandedMatrix<double> & AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M8>::value(BandedMatrix<double> & m2, BandedMatrix<double> & result, DenseVector<double> & d, unsigned long d_width, unsigned long d_height)
{
    DenseVector<double> m8diag(m2.size());
    DenseVector<double> m8bandplus3(m2.size());
    DenseVector<double> m8bandplus6(m2.size());
    DenseVector<double> m8bandminus3(m2.size());

    double * m8diag_e = m8diag.elements();
    double * m8bandplus3_e = m8bandplus3.elements();
    double * m8bandplus6_e = m8bandplus6.elements();
    double * m8bandminus3_e = m8bandminus3.elements();

    double * m2diag_e = (m2.band(0)).elements();
    double * m2bandplus3_e = (m2.band(3)).elements();
    double * m2bandplus6_e = (m2.band(6)).elements();
    double * m2bandminus3_e = (m2.band(-3)).elements();

    unsigned long tripels(m2.size()/3);
    unsigned long sse_tripels((tripels - (tripels % 2)) / 2);
    unsigned long sse_limit(sse_tripels * 6);

    double HONEI_ATTRIBUTE(aligned(16)) d1_buffer[2];
    double HONEI_ATTRIBUTE(aligned(16)) d2_buffer[2];
    double HONEI_ATTRIBUTE(aligned(16)) d3_buffer[2];

    for(int i(0); i < 2; ++i)
    {
        d1_buffer[i] = d[0];
        d2_buffer[i] = d[1];
        d3_buffer[i] = d[2];
    }

    double HONEI_ATTRIBUTE(aligned(16)) d_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) d_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) d_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) b1_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b1_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b1_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) b2_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b2_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) b2_q2_buffer[tripels];

    double HONEI_ATTRIBUTE(aligned(16)) bm1_h_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) bm1_q1_buffer[tripels];
    double HONEI_ATTRIBUTE(aligned(16)) bm1_q2_buffer[tripels];

    ///Alignment and compaction loop
    unsigned long buffer_i(0);
    for(unsigned long i(0); i < m2.size(); i += 3)
    {
        d_h_buffer[buffer_i] = m2diag_e[i];
        d_q1_buffer[buffer_i] = m2diag_e[i + 1];
        d_q2_buffer[buffer_i] = m2diag_e[i + 2];
        b1_h_buffer[buffer_i] = m2bandplus3_e[i];
        b1_q1_buffer[buffer_i] = m2bandplus3_e[i + 1];
        b1_q2_buffer[buffer_i] = m2bandplus3_e[i + 2];
        b2_h_buffer[buffer_i] = m2bandplus6_e[i];
        b2_q1_buffer[buffer_i] = m2bandplus6_e[i + 1];
        b2_q2_buffer[buffer_i] = m2bandplus6_e[i + 2];
        bm1_h_buffer[buffer_i] = m2bandminus3_e[i];
        bm1_q1_buffer[buffer_i] = m2bandminus3_e[i + 1];
        bm1_q2_buffer[buffer_i] = m2bandminus3_e[i + 2];
        ++buffer_i;
    }

    __m128d m0, m1, m3, m4, m5, m6, m7, m8;

    m0 = _mm_load_pd(d1_buffer);
    m1 = _mm_load_pd(d2_buffer);
    m3 = _mm_load_pd(d3_buffer);

    m0 = _mm_mul_pd(m0, m0);
    m1 = _mm_mul_pd(m1, m1);
    m3 = _mm_mul_pd(m3, m3);

    _mm_store_pd(d1_buffer, m0);
    _mm_store_pd(d2_buffer, m1);
    _mm_store_pd(d3_buffer, m3);


    for(unsigned long i(0); i < tripels - (tripels % 2); i += 2)
    {
        m4 = _mm_load_pd(d_h_buffer + i);
        m5 = _mm_load_pd(b1_h_buffer + i);
        m6 = _mm_load_pd(b2_h_buffer + i);
        m7 = _mm_load_pd(bm1_h_buffer + i);

        m4 = _mm_mul_pd(m4, m0);
        m5 = _mm_mul_pd(m5, m0);
        m6 = _mm_mul_pd(m6, m0);
        m7 = _mm_mul_pd(m7, m0);

        _mm_store_pd(d_h_buffer + i, m4);
        _mm_store_pd(b1_h_buffer + i, m5);
        _mm_store_pd(b2_h_buffer + i, m6);
        _mm_store_pd(bm1_h_buffer + i, m7);

        m4 = _mm_load_pd(d_q1_buffer + i);
        m5 = _mm_load_pd(b1_q1_buffer + i);
        m6 = _mm_load_pd(b2_q1_buffer + i);
        m7 = _mm_load_pd(bm1_q1_buffer + i);

        m4 = _mm_mul_pd(m4, m1);
        m5 = _mm_mul_pd(m5, m1);
        m6 = _mm_mul_pd(m6, m1);
        m7 = _mm_mul_pd(m7, m1);

        _mm_store_pd(d_q1_buffer + i, m4);
        _mm_store_pd(b1_q1_buffer + i, m5);
        _mm_store_pd(b2_q1_buffer + i, m6);
        _mm_store_pd(bm1_q1_buffer + i, m7);

        m4 = _mm_load_pd(d_q2_buffer + i);
        m5 = _mm_load_pd(b1_q2_buffer + i);
        m6 = _mm_load_pd(b2_q2_buffer + i);
        m7 = _mm_load_pd(bm1_q2_buffer + i);

        m4 = _mm_mul_pd(m4, m3);
        m5 = _mm_mul_pd(m5, m3);
        m6 = _mm_mul_pd(m6, m3);
        m7 = _mm_mul_pd(m7, m3);

        _mm_store_pd(d_q2_buffer + i, m4);
        _mm_store_pd(b1_q2_buffer + i, m5);
        _mm_store_pd(b2_q2_buffer + i, m6);
        _mm_store_pd(bm1_q2_buffer + i, m7);
    }
    ///Result loop for SSE computed elements
    buffer_i = 0;
    for(unsigned long i(0); i < sse_limit; i += 3)
    {
        m8diag_e[i] = d_h_buffer[buffer_i];
        m8diag_e[i + 1] = d_q1_buffer[buffer_i];
        m8diag_e[i + 2] = d_q2_buffer[buffer_i];
        m8bandplus3_e[i] = b1_h_buffer[buffer_i];
        m8bandplus3_e[i + 1] = b1_q1_buffer[buffer_i];
        m8bandplus3_e[i + 2] = b1_q2_buffer[buffer_i];
        m8bandplus6_e[i] = b2_h_buffer[buffer_i];
        m8bandplus6_e[i + 1] = b2_q1_buffer[buffer_i];
        m8bandplus6_e[i + 2] = b2_q2_buffer[buffer_i];
        m8bandminus3_e[i] = bm1_h_buffer[buffer_i];
        m8bandminus3_e[i + 1] = bm1_q1_buffer[buffer_i];
        m8bandminus3_e[i + 2] = bm1_q2_buffer[buffer_i];
        ++buffer_i;
    }

    ///Non SSE loop
    for(unsigned long i(sse_limit); i < m2.size(); i += 3)
    {
        m8diag_e[i] = d1_buffer[0] * m2diag_e[i];
        m8diag_e[i + 1] = d2_buffer[0] * m2diag_e[i + 1];
        m8diag_e[i + 2] = d3_buffer[0] * m2diag_e[i + 2];
        m8bandplus3_e[i] = d1_buffer[0] * m2bandplus3_e[i];
        m8bandplus3_e[i + 1] = d2_buffer[0] * m2bandplus3_e[i + 1];
        m8bandplus3_e[i + 2] = d3_buffer[0] * m2bandplus3_e[i + 2];
        m8bandplus6_e[i] = d1_buffer[0] * m2bandplus6_e[i];
        m8bandplus6_e[i + 1] = d2_buffer[0] * m2bandplus6_e[i + 1];
        m8bandplus6_e[i + 2] = d3_buffer[0] * m2bandplus6_e[i + 2];
        m8bandminus3_e[i] = d1_buffer[0] * m2bandminus3_e[i];
        m8bandminus3_e[i + 1] = d2_buffer[0] * m2bandminus3_e[i + 1];
        m8bandminus3_e[i + 2] = d3_buffer[0] * m2bandminus3_e[i + 2];
    }

    result.insert_band(0, m8diag);
    result.insert_band(3, m8bandplus3);
    result.insert_band(6, m8bandplus6);
    result.insert_band(-3, m8bandminus3);

    return result;
}

