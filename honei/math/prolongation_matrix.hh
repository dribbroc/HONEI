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

#ifndef MATH_GUARD_PROLONGATION_MATRIX_HH
#define MATH_GUARD_PROLONGATION_MATRIX_HH 1

#include <honei/la/sparse_matrix.hh>
#include <cmath>
#include <iostream>

namespace honei
{
    ///for squared regular grids only
    template <typename Tag_>
    class ProlongationMatrix
    {
        private:

            static inline unsigned long _coarse_index(unsigned long fine_index, unsigned long fine_grid_row, unsigned long offset)
            {
                return (fine_index / 2) - ((fine_grid_row / 2) * ((offset - 1) / 2));
            }

        public:

            template <typename DT_>
            static void value(SparseMatrix<DT_> & target, DenseVector<unsigned long> & indices_fine, DenseVector<unsigned long> & indices_coarse)
            {
                CONTEXT("During assembly of prolongation matrix: ");

                if(target.rows() != indices_fine.size())
                    throw InternalError("Matrix row count " + stringify(target.rows()) + " does not match size of fine index array (" + stringify(indices_fine.size()) + ")");
                if(target.columns() != indices_coarse.size())
                    throw InternalError("Matrix column count " + stringify(target.columns()) + " does not match size of coarse index array (" + stringify(indices_coarse.size()) + ")");
                unsigned long fine_grid_dim((unsigned long)sqrt(indices_fine.size())); ///we have a squared grid

                for(unsigned long i(0) ; i < fine_grid_dim ; ++i)
                    for(unsigned long j(0) ; j < fine_grid_dim ; ++j)
                    {
                        DT_ coeff(0);
                        unsigned long i_fine, i_coarse;
                        i_fine = i * fine_grid_dim + j;

                        if(i % 2 == 0)
                        {
                            if(j % 2 == 0) ///exact on coarse grid -> write only one coeff
                            {
                                coeff = DT_(1);
                                i_coarse = _coarse_index(i_fine, i, fine_grid_dim);
                                target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            }
                            else ///on coarse grid row but not on coarse grid column -> write two coeffs
                            {
                                coeff = DT_(0.5);
                                i_coarse = _coarse_index(i_fine - 1, i, fine_grid_dim);
                                target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                                i_coarse = _coarse_index(i_fine + 1, i, fine_grid_dim);
                                target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            }
                        }
                        else if(j % 2 == 0) ///on coarse grid column but not on row -> write two coeffs
                        {
                            coeff = DT_(0.5);
                            i_coarse = _coarse_index(i_fine + fine_grid_dim, i + 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            i_coarse = _coarse_index(i_fine - fine_grid_dim, i - 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                        }
                        else ///neither on coarse grid row nor column -> write four coeffs
                        {
                            coeff = DT_(0.25);
                            i_coarse = _coarse_index(i_fine + fine_grid_dim + 1, i + 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            i_coarse = _coarse_index(i_fine + fine_grid_dim - 1, i + 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            i_coarse = _coarse_index(i_fine - fine_grid_dim + 1, i - 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                            i_coarse = _coarse_index(i_fine - fine_grid_dim - 1, i - 1, fine_grid_dim);
                            target(indices_fine[i_fine] , indices_coarse[i_coarse]) = coeff;
                        }
                    }
            }
    };
}

#endif
