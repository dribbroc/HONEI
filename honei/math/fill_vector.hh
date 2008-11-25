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

#ifndef MATH_GUARD_FILL_VECTOR_HH
#define MATH_GUARD_FILL_VECTOR_HH 1


#include<honei/la/dense_vector.hh>
#include<honei/la/banded_matrix_q1.hh>
#include<honei/math/methods.hh>
#include<cmath>
#include <iostream>

using namespace methods;
namespace honei
{
    template <typename Tag_, typename Application_, typename BoundaryType_>
        struct FillVector
        {
        };

    template <typename Tag_>
        struct FillVector<Tag_, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>
        {
            private:
                template <typename DT_>
                    static inline DT_ _f_1(DT_ y)
                    {
                        //f_1(y) = 1 + 2y, y in [0,1/2)
                        //       = 2,      y = 0
                        //       = 3  - 2y,y in (1/2, 1]

                        if(y < 0.5)
                            return DT_(1. + 2 * y);
                        else if(y > 0.5)
                            return DT_(3. - 2.* y);
                        else
                            return DT_(2.);
                    }
            public:
                template <typename DT_>
                    static inline void value(DenseVector<DT_> & target)
                    {
                        unsigned long size(target.size());
                        unsigned long N(size);
                        unsigned long M((unsigned long)sqrt(N));
                        /*for (unsigned long i(1) ; i <= N ; ++i)
                        {
                            // first, Dirichlet unit rows
                            if (i <= M || i > N-M || (i%M) == 1)
                            {
                                target[i-1] = 0.0;
                            }
                            // then, Neumann on the right
                            else if ((i % M) == 0 && i > M && i < N)
                            {
                                target[i-1] = 1.;//_f_1((i - M)/M * DT_(1./(DT_(M) - 1.)));
                            }
                            // then, inner points
                            else
                            {
                                target[i-1] = 1.;
                            }
                        }*/

                        unsigned long level(0);

                        switch(N)
                        {
                            case 1050625:
                                {
                                    level = 10;
                                }
                                break;
                            case 263169:
                                {
                                    level = 9;
                                }
                                break;
                            case 66049:
                                {
                                    level = 8;
                                }
                                break;
                            case 16641:
                                {
                                    level = 7;
                                }
                                break;
                            case 4225:
                                {
                                    level = 6;
                                }
                                break;
                            case 1089:
                                {
                                    level = 5;
                                }
                                break;
                            case 289:
                                {
                                    level = 4;
                                }
                                break;
                            case 81:
                                {
                                    level = 3;
                                }
                                break;
                            case 25:
                                {
                                    level = 2;
                                }
                                break;
                            case 9:
                                {
                                    level = 1;
                                }
                                break;
                        }
                        for (unsigned long i(1); i <= N; i++)
                        {
                            // first, Dirichlet
                            if (i <= M || i > N - M || (i%M) == 1)
                            {
                                target[i-1] = 2.0;
                            }
                            // then, Neumann
                            else if ((i%M) == 0 && i>M && i<N)
                            {
                                target[i-1] = 0.5*(1.0/pow(4.0, level));
                            }
                            // then, inner points
                            else
                            {
                                target[i-1] = 1.0/pow(4.0, level);
                            }
                        }
                    }
        };
}

#endif
