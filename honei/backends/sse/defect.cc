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
                unsigned long stride, unsigned long rows, unsigned long /*num_cols_per_row*/, const unsigned long threads)
        {
            for (unsigned long row(0) ; row < rows ; ++row)
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

                        //if (A_ij != 0)
                        {
                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }
                    }

                    tAj += stride;
                    tAx += stride;
                }
                result[row] = rhs[row] - sum;
            }
        }

        void defect_smell_dv(float * result, const float * rhs, const unsigned long * Aj, const float * Ax, const unsigned long * Arl, const float * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
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

                        //if (A_ij != 0)
                        {
                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }
                    }

                    tAj += stride;
                    tAx += stride;
                }
                result[row] = rhs[row] - sum;
            }
        }

        void defect_smell_dv(double * result, const double * rhs, const unsigned long * Aj, const double * Ax, const unsigned long * Arl, const double * b,
                unsigned long stride, unsigned long rows, unsigned long /*num_cols_per_row*/, const unsigned long threads)
        {
            for (unsigned long row(0) ; row < rows ; ++row)
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

                        //if (A_ij != 0)
                        {
                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }
                    }

                    tAj += stride;
                    tAx += stride;
                }
                result[row] = rhs[row] - sum;
            }
        }

        void defect_smell_dv(double * result, const double * rhs, const unsigned long * Aj, const double * Ax, const unsigned long * Arl, const double * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
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

                        //if (A_ij != 0)
                        {
                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }
                    }

                    tAj += stride;
                    tAx += stride;
                }
                result[row] = rhs[row] - sum;
            }
        }
    }
}
