/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008-2010 Markus Geveler <apryde@gmx.de>
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

#pragma once
#ifndef MATH_GUARD_RESTRICT_HH
#define MATH_GUARD_RESTRICT_HH 1

#include<honei/la/dense_vector.hh>
#include<honei/la/banded_matrix.hh>
#include<honei/la/sparse_matrix.hh>
#include<honei/la/product.hh>
#include<honei/util/attributes.hh>
#include<honei/math/methods.hh>
#include<honei/math/apply_dirichlet_boundaries.hh>
#include<cmath>

namespace honei
{
    template<typename Tag_, typename Type_>
    struct Restriction
    {
    };

    template<typename Tag_>
    struct Restriction<Tag_, methods::NONE>
    {
        public:
            template <typename Prec_, typename MatrixType_>
                static DenseVector<Prec_> & value(DenseVector<Prec_>&  coarse, DenseVector<Prec_>& fine, DenseVector<unsigned long>& mask, HONEI_UNUSED MatrixType_ & resmat)
                {
                    unsigned long n_fine(fine.size());
                    unsigned long n_x_fine((unsigned long)sqrt((Prec_)n_fine) - 1);
                    unsigned long n_x_coarse(n_x_fine / 2);
                    unsigned long n_y_fine(n_x_fine + 1);
                    unsigned long n_y_coarse(n_x_coarse + 1);
                    unsigned long i_step(2*n_y_fine);

                    // temps
                    int i_aux, i_aux_1;

                    // x___________x
                    // |     |     |
                    // |     |     |
                    // |_____|_____|
                    // |     |     |
                    // |     |     |
                    // x_____|_____x
                    for (unsigned long i(1) ; i <= n_y_coarse ; ++i)
                    {
                        for(unsigned long j(0), k(0) ; j < n_y_coarse ; ++j, k +=2)
                        {
                            coarse[1 + n_y_coarse * (i - 1) - 1 + j] = fine[2 * n_y_fine * (i - 1) + k];
                        }
                    }

                    for(unsigned i(1); i <= n_x_coarse; ++i)
                    {
                        i_aux = 1 + n_y_fine + 2 * n_y_fine * (i - 1);
                        i_aux_1 = 1 + n_y_coarse * i;

                        // x_______________________x
                        // | _         |         _ |
                        // ||\         |         /||
                        // | |         |         | |
                        // | |1/2      |     1/2 | |
                        // | |         |         | |
                        // | /         |         \ |
                        // o___________|___________o
                        // |           |           |
                        // | \         |         / |
                        // | |         |         | |
                        // | |1/2      |     1/2 | |
                        // | |         |         | |
                        // ||/_        |        _\||
                        // x___________|___________x

                        for(unsigned long j(0) ; j < n_y_coarse; ++j)
                        {
                            coarse[i_aux_1 - 1 + j] += 0.5 * fine[i_aux - 1 + j * 2];
                        }

                        i_aux_1 = i_aux_1 - n_y_coarse;

                        for(unsigned long j(0) ; j < n_y_coarse; ++j)
                        {
                            coarse[i_aux_1 - 1 + j] += 0.5 * fine[i_aux - 1 + j * 2];
                        }

                        // x___________o___________x
                        // | \_______/ | \_______/ |
                        // |    1/2    |    1/2    |
                        // |           |           |
                        // |           |           |
                        // |           |           |
                        // |           |           |
                        // |___________|___________|
                        // |           |           |
                        // |           |           |
                        // |           |           |
                        // |           |           |
                        // |  _______  |  _______  |
                        // | /  1/2  \ | /  1/2  \ |
                        // x___________o___________x

                        for(unsigned long j(0) ; j < n_y_coarse; ++j)
                        {
                            coarse[i - 1 + j * n_y_coarse] += 0.5 * fine[2 * i - 1 + j * i_step];
                        }

                        for(unsigned long j(0) ; j < n_y_coarse; ++j)
                        {
                            coarse[i + j * n_y_coarse] += 0.5 * fine[2 * i - 1 + j * i_step];
                        }

                        ++i_aux;

                        // _______________________
                        // |  _        |         _ |
                        // | |\        |         /||
                        // |   \ 1/4   |    1/4/   |
                        // |     \     |     /     |
                        // |       \   |   /       |
                        // |         \ | /         |
                        // |___________o___________|
                        // |         / | \         |
                        // |        /  |  \        |
                        // |      /1/4 |  1/4\     |
                        // |     /     |      \    |
                        // |   /       |        \  |
                        // | |/_       |        _\||
                        // x___________|___________x

                        for(unsigned long j(0) ; j < n_x_coarse; ++j)
                        {
                            coarse[i_aux_1 -1 + j] += 0.25 * fine[i_aux - 1 + j * 2];
                        }

                        for(unsigned long j(0) ; j < n_x_coarse; ++j)
                        {
                            coarse[i_aux_1 + j] += 0.25 * fine[i_aux - 1 + j * 2];
                        }

                        for(unsigned long j(0) ; j < n_x_coarse; ++j)
                        {
                            coarse[n_y_coarse * i + j] += 0.25 * fine[i_aux - 1 + j * 2];
                        }

                        for(unsigned long j(0) ; j < n_x_coarse; ++j)
                        {
                            coarse[1 + n_y_coarse * i + j] += 0.25 * fine[i_aux - 1 + j * 2];
                        }
                    }
                    ApplyDirichletBoundaries<Tag_>::value(coarse, mask);
                    return coarse;
                }
    };

    template<typename Tag_>
    struct Restriction<Tag_, methods::PROLMAT>
    {
        public:
            template <typename Prec_, typename MatrixType_>
                static DenseVector<Prec_> & value(DenseVector<Prec_>&  coarse, DenseVector<Prec_>& fine, HONEI_UNUSED DenseVector<unsigned long>& mask, HONEI_UNUSED MatrixType_ & resmat)
                {
                    Product<Tag_>::value(coarse, resmat, fine);
                    //ApplyDirichletBoundaries<Tag_>::value(coarse, mask);
                    return coarse;
                }
    };

    template <> struct Restriction<tags::GPU::CUDA, methods::NONE>
    {
        static DenseVector<float> & value(DenseVector<float> & coarse,
                const DenseVector<float> & fine, const DenseVector<unsigned long> & mask, HONEI_UNUSED BandedMatrixQ1<float> & resmat);

        static DenseVector<double> & value(DenseVector<double> & coarse,
                const DenseVector<double> & fine, const DenseVector<unsigned long> & mask, HONEI_UNUSED BandedMatrixQ1<double> & resmat);

        static DenseVector<float> & value(DenseVector<float> & coarse,
                const DenseVector<float> & fine, const DenseVector<unsigned long> & mask, HONEI_UNUSED SparseMatrixELL<float> & resmat);

        static DenseVector<double> & value(DenseVector<double> & coarse,
                const DenseVector<double> & fine, const DenseVector<unsigned long> & mask, HONEI_UNUSED SparseMatrixELL<double> & resmat);
    };
}
#endif
