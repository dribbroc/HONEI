/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LBM_GUARD_PARTIAL_DERIVATIVE_HH
#define LBM_GUARD_PARTIAL_DERIVATIVE_HH 1

/**
 * \file
 * Implementation of partial derivative modules used by  LBM - (SWE) solvers using a
 * PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <cmath>

using namespace honei::lbm;
using namespace honei::lbm::lbm_directions;
using namespace honei::lbm::lbm_source_schemes;

namespace honei
{
    template<typename Tag_, typename Direction_, typename Discretization_>
        class PartialDerivative
        {
        };

    template<typename Tag_>
        class PartialDerivative<Tag_, X, CENTRALDIFF>
        {
            public:
                template<typename DT_>
                    static DenseMatrix<DT_> value(DenseMatrix<DT_>& f, DT_ delta_x)
                    {
                        DenseMatrix<DT_> result(f.rows(), f.columns());

                        for(unsigned long i(0) ; i < f.rows() ; ++i)
                        {
                            for(unsigned long j(0) ; j < f.columns() ; ++j)
                            {
                                unsigned long left, right;
                                if(j == 0)
                                    left = j;
                                else
                                    left = j - 1;
                                if(j == f.columns() - 1)
                                    right = j;
                                else
                                    right = j + 1;

                                result(i , j) = (f(i , right) - f(i , left))/(DT_(2) * delta_x);
                            }
                        }
                        return result;
                    }
        };

    template<typename Tag_>
        class PartialDerivative<Tag_, Y, CENTRALDIFF>
        {
            public:
                template<typename DT_>
                    static DenseMatrix<DT_> value(DenseMatrix<DT_>& f, DT_ delta_y)
                    {
                        DenseMatrix<DT_> result(f.rows(), f.columns());

                        for(unsigned long i(0) ; i < f.rows() ; ++i)
                        {
                            for(unsigned long j(0) ; j < f.columns() ; ++j)
                            {
                                unsigned long upper, lower;
                                if(i == 0)
                                    upper = i;
                                else
                                    upper = i - 1;
                                if(i == f.columns() - 1)
                                    lower = i;
                                else
                                    lower = i + 1;

                                result(i , j) = (f(upper , j) - f(lower , j))/(DT_(2) * delta_y);
                            }
                        }
                        return result;
                    }
        };
}
#endif
