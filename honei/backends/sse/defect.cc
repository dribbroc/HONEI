/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/attributes.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace sse
    {
        void defect_smell_dv(float * result, const float * rhs, const unsigned long * Aj, const float * Ax, const unsigned long * Arl, const float * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
        {
            if (threads % 4 != 0)
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const float * tAx(Ax);
                    float sum(0);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; ++thread)
                        {
                            const float A_ij = *(tAx + thread);

                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    result[row] = rhs[row] - sum;
                }
            }
            else
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const float * tAx(Ax);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    __m128 A_ij, b_v, p;
                    union sse4
                    {
                        __m128 m;
                        float f[4];
                    } sum_v;
                    sum_v.m = _mm_setzero_ps();

                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; thread+=4)
                        {
                            A_ij = _mm_load_ps(tAx + thread);
                            const unsigned long col0 = *(tAj + thread);
                            const unsigned long col1 = *(tAj + thread + 1);
                            const unsigned long col2 = *(tAj + thread + 2);
                            const unsigned long col3 = *(tAj + thread + 3);
                            b_v = _mm_set_ps(*(b+col3), *(b+col2), *(b+col1), *(b+col0));
                            p = _mm_mul_ps(A_ij, b_v);
                            sum_v.m = _mm_add_ps(p, sum_v.m);
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    float sum = sum_v.f[0] + sum_v.f[1] + sum_v.f[2] + sum_v.f[3];

                    result[row] = rhs[row] - sum;
                }
            }
        }

        void defect_smell_dv(double * result, const double * rhs, const unsigned long * Aj, const double * Ax, const unsigned long * Arl, const double * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
        {
            if (threads % 2 != 0)
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const double * tAx(Ax);
                    double sum(0);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; ++thread)
                        {
                            const double A_ij = *(tAx + thread);

                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    result[row] = rhs[row] - sum;
                }
            }
            else
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const double * tAx(Ax);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    __m128d A_ij, b_v, p;
                    union sse2
                    {
                        __m128d m;
                        double d[2];
                    } sum_v;
                    sum_v.m = _mm_setzero_pd();

                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; thread+=2)
                        {
                            A_ij = _mm_load_pd(tAx + thread);
                            const unsigned long col0 = *(tAj + thread);
                            const unsigned long col1 = *(tAj + thread + 1);
                            b_v = _mm_set_pd(*(b+col1), *(b+col0));
                            p = _mm_mul_pd(A_ij, b_v);
                            sum_v.m = _mm_add_pd(p, sum_v.m);
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    double sum = sum_v.d[0] + sum_v.d[1];

                    result[row] = rhs[row] - sum;
                }
            }
        }

        void defect_csr_dv(float * result, const float * rhs, const unsigned long * Aj, const float * Ax, const unsigned long * Ar, const float * b,
                unsigned long blocksize, unsigned long row_start, unsigned long row_end)
        {
            if (blocksize != 4)
            {
                switch(blocksize)
                {
                    case 1:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                sum += Ax[i] * b[Aj[i]];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    case 2:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    /*case 4:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[4*i] * b[col];
                                sum += Ax[4*i+1] * b[col+1];
                                sum += Ax[4*i+2] * b[col+2];
                                sum += Ax[4*i+3] * b[col+3];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;*/
                    case 8:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[8*i] * b[col];
                                sum += Ax[8*i+1] * b[col+1];
                                sum += Ax[8*i+2] * b[col+2];
                                sum += Ax[8*i+3] * b[col+3];
                                sum += Ax[8*i+4] * b[col+4];
                                sum += Ax[8*i+5] * b[col+5];
                                sum += Ax[8*i+6] * b[col+6];
                                sum += Ax[8*i+7] * b[col+7];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    default:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                                {
                                    sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                                }
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                }
            }
            else
            {
                __m128 A_ij, b_v;
                union sse4
                {
                    __m128 m;
                    float f[4];
                } sum_v;

                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    sum_v.m = _mm_setzero_ps();

                    const unsigned long end(Ar[row+1]);
                    for (unsigned long i(Ar[row]) ; i < end ; i++)
                    {
                        A_ij = _mm_load_ps(Ax + (i * blocksize));
                        const unsigned long col = Aj[i];
                        b_v = _mm_load_ps(b+col);
                        b_v = _mm_mul_ps(A_ij, b_v);
                        sum_v.m = _mm_add_ps(b_v, sum_v.m);
                    }
                    float sum = sum_v.f[0] + sum_v.f[1] + sum_v.f[2] + sum_v.f[3];

                    result[row] = rhs[row] - sum;
                }
            }
        }

        void defect_csr_dv(double * result, const double * rhs, const unsigned long * Aj, const double * Ax, const unsigned long * Ar, const double * b,
                unsigned long blocksize, unsigned long row_start, unsigned long row_end)
        {
            if (blocksize != 2)
            {
                switch(blocksize)
                {
                    case 1:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                sum += Ax[i] * b[Aj[i]];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    /*case 2:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;*/
                    case 4:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[4*i] * b[col];
                                sum += Ax[4*i+1] * b[col+1];
                                sum += Ax[4*i+2] * b[col+2];
                                sum += Ax[4*i+3] * b[col+3];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    case 8:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[8*i] * b[col];
                                sum += Ax[8*i+1] * b[col+1];
                                sum += Ax[8*i+2] * b[col+2];
                                sum += Ax[8*i+3] * b[col+3];
                                sum += Ax[8*i+4] * b[col+4];
                                sum += Ax[8*i+5] * b[col+5];
                                sum += Ax[8*i+6] * b[col+6];
                                sum += Ax[8*i+7] * b[col+7];
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                    default:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                                {
                                    sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                                }
                            }
                            result[row] = rhs[row] - sum;
                        }
                        break;
                }
            }
            else
            {
                __m128d A_ij, b_v;
                union sse4
                {
                    __m128d m;
                    double d[4];
                } sum_v;

                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    sum_v.m = _mm_setzero_pd();

                    const unsigned long end(Ar[row+1]);
                    for (unsigned long i(Ar[row]) ; i < end ; i++)
                    {
                        A_ij = _mm_load_pd(Ax + (i * blocksize));
                        const unsigned long col = Aj[i];
                        b_v = _mm_load_pd(b+col);
                        b_v = _mm_mul_pd(A_ij, b_v);
                        sum_v.m = _mm_add_pd(b_v, sum_v.m);
                    }
                    double sum = sum_v.d[0] + sum_v.d[1];

                    result[row] = rhs[row] - sum;
                }
            }
        }
    }
}
