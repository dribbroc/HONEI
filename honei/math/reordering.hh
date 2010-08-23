/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI math C++ library. HONEI is free software;
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

#pragma once
#ifndef MATH_GUARD_REORDERING_HH
#define MATH_GUARD_REORDERING_HH 1

#include<honei/math/methods.hh>
#include<honei/la/dense_vector.hh>
#include<cmath>

namespace honei
{
    template <typename Tag_, typename OrderType_>
    class Reordering
    {
    };

    template <typename Tag_>
    class Reordering<Tag_, methods::NATURAL>
    {
        public:

        static void value(DenseVector<unsigned long> & fine, DenseVector<unsigned long> & coarse)
        {
            for(unsigned long i(0) ; i < fine.size() ; ++i)
                fine[i] = i;

            for(unsigned long i(0) ; i < coarse.size() ; ++i)
                coarse[i] = i;
        }
    };

    template <typename Tag_>
    class Reordering<Tag_, methods::TWO_LEVEL>
    {
        public:

        static void value(DenseVector<unsigned long> & fine, DenseVector<unsigned long> & coarse)
        {
            unsigned long fine_grid_dim((unsigned long)sqrt(fine.size()));
            unsigned long fine_index(coarse.size());
            unsigned long coarse_index(0);
            for(unsigned long i(0) ; i < fine_grid_dim ; ++i)
                for(unsigned long j(0) ; j < fine_grid_dim ; ++j)
                {
                    if(i % 2 == 0 && j % 2 == 0)
                    {
                        fine[i * fine_grid_dim + j] = coarse_index;
                        ++coarse_index;
                    }
                    else
                    {
                        fine[i * fine_grid_dim + j] = fine_index;
                        ++fine_index;
                    }
                }

            for(unsigned long i(0) ; i < coarse.size() ; ++i)
                coarse[i] = i;
        }
    };
}

#endif
