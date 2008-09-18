/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MATH_GUARD_APPLY_DIRICHLET_BOUNDARIES_HH
#define MATH_GUARD_APPLY_DIRICHLET_BOUNDARIES_HH 1

#include <honei/la/dense_vector.hh>
#include <cmath>

namespace honei
{
    template<typename Tag_>
        class ApplyDirichletBoundaries
        {
            public:
                template<typename Prec_>
                    static DenseVector<Prec_>& value(DenseVector<Prec_>& grid, DenseVector<unsigned long>& mask)
                    {
                        int i;
                        int N = grid.size();
                        int M = (int)sqrt((double)N);

                        // bottom left node
                        if (mask[0] == 2)
                            grid[0] = 0.0f;
                        // bottom right node
                        if (mask[1] == 2)
                            grid[M-1] = 0.0f;
                        // top right node
                        if (mask[2] == 2)
                            grid[N-1] = 0.0f;
                        // top left node
                        if (mask[3] == 2)
                            grid[N-M] = 0.0f;
                        // bottom edge
                        if (mask[4] == 2)
                        {
                            for (unsigned long i(1) ; i < M - 1 ; i++)
                            {
                                grid[i] = 0.0f;
                            }
                        }
                        // right edge
                        if (mask[5] == 2)
                        {
                            for (unsigned long i(2 * M - 1) ; i < N - M ; i += M)
                            {
                                grid[i] = 0.0f;
                            }
                        }
                        // top edge
                        if (mask[6] == 2)
                        {
                            for (unsigned long i(N - M + 1) ; i < N - 1; i++)
                            {
                                grid[i] = 0.0f;
                            }
                        }
                        // left edge
                        if (mask[7] == 2)
                        {
                            for (unsigned long i(M) ; i < N - M; i += M)
                            {
                                grid[i] = 0.0f;
                            }
                        }
                    }
        };
}

#endif
