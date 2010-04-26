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

#ifndef MATH_GUARD_PROLONGATE_HH
#define MATH_GUARD_PROLONGATE_HH 1

#include<honei/la/dense_vector.hh>
#include<honei/la/banded_matrix.hh>
#include<honei/la/sparse_matrix.hh>
#include<honei/la/product.hh>
#include<honei/util/attributes.hh>
#include<honei/math/methods.hh>
#include<honei/math/apply_dirichlet_boundaries.hh>
#include<cmath>

using namespace methods;
namespace honei
{
    template<typename Tag_, typename Type_>
        class Prolongation
        {
        };

    template<typename Tag_>
        class Prolongation<Tag_, NONE>
        {
            public:
                template <typename Prec_, typename MatrixType_>
                    static DenseVector<Prec_> & value(DenseVector<Prec_>&  fine,
                            DenseVector<Prec_>& coarse,
                            DenseVector<unsigned long>& mask,
                            HONEI_UNUSED MatrixType_ & prolmat)
                    {
                        unsigned long n_fine(fine.size());
                        unsigned long n_x_fine((unsigned long)sqrt((Prec_)n_fine) - 1);
                        unsigned long n_x_coarse(n_x_fine / 2);
                        unsigned long n_y_fine(n_x_fine + 1);
                        unsigned long n_y_coarse(n_x_coarse + 1);
                        unsigned long i_step(2*n_y_fine);

                        for (unsigned long i(0); i < n_fine; ++i)
                        {
                            fine[i] = 0.0;
                        }

                        for (unsigned long i(1); i <= n_y_coarse; ++i)
                        {
                            // x___________x
                            // |     |     |
                            // |     |     |
                            // |_____|_____|
                            // |     |     |
                            // |     |     |
                            // x_____|_____x

                            for (unsigned long j(0), k(0); j < n_y_coarse; ++j, k+=2)
                            {
                                fine[i_step*(i-1) + k] = coarse[i + n_x_coarse * (i-1) - 1 + j];
                            }
                        }

                        for (unsigned long  i(1); i <= n_x_coarse; ++i)
                        {
                            //  ___________
                            // |     |     |
                            // |     |     |
                            // x_____|_____x
                            // |     |     |
                            // |     |     |
                            // |_____|_____|

                            // x_new = 1/2 * (x_b + x_t)
                            unsigned long iaux(2 * i);

                            for(unsigned long j(0); j < n_y_coarse; ++j)
                            {
                                fine[iaux - 1 + j*i_step] += coarse[i - 1 + j * n_y_coarse] * 0.5;
                                fine[iaux - 1 + j*i_step] += coarse[i + j * n_y_coarse] * 0.5;
                            }
                            //  _____x_____
                            // |     |     |
                            // |     |     |
                            // |_____|_____|
                            // |     |     |
                            // |     |     |
                            // |_____x_____|

                            // x_new = 1/2 * (x_l + x_r)
                            iaux = n_y_fine + 1 + (i - 1) * i_step;

                            for(unsigned long j(0); j < n_y_coarse; ++j)
                            {
                                fine[iaux - 1 + j * 2] += fine[iaux - n_y_fine - 1 + j * 2] * 0.5;
                                fine[iaux - 1 + j * 2] += fine[iaux + n_y_fine - 1+ j * 2] * 0.5;
                            }
                        }

                        for (unsigned long i(1); i <= n_x_coarse; ++i)
                        {
                            //  ___________
                            // |     |     |
                            // |     |     |
                            // |_____x_____|
                            // |     |     |
                            // |     |     |
                            // |_____|_____|

                            unsigned long iaux(n_x_fine + 3 + (i - 1) * i_step);

                            // x_neq = 1/4 * (x_lb + x_rb + x_lt + x_rt))
                            for(unsigned long j(0) ; j < n_x_coarse ; ++j)
                            {
                                fine[iaux - 1 + j * 2] += fine[iaux - n_x_fine - 3 + j * 2] * 0.25;
                                fine[iaux - 1 + j * 2] += fine[iaux - n_x_fine - 1 + j * 2] * 0.25;
                                fine[iaux - 1 + j * 2] += fine[iaux + n_x_fine - 1 + j * 2] * 0.25;
                                fine[iaux - 1 + j * 2] += fine[iaux + n_x_fine + 1 + j * 2] * 0.25;
                            }
                        }
                        ApplyDirichletBoundaries<Tag_>::value(fine, mask);
                        return fine;
                    }
        };

    template<typename Tag_>
        class Prolongation<Tag_, PROLMAT>
        {
            public:
                template <typename Prec_, typename MatrixType_>
                    static DenseVector<Prec_> & value(DenseVector<Prec_>&  fine,
                            DenseVector<Prec_>& coarse,
                            DenseVector<unsigned long>& mask,
                            HONEI_UNUSED MatrixType_ & prolmat)
                    {
                        Product<Tag_>::value(fine, prolmat, coarse);
                        //ApplyDirichletBoundaries<Tag_>::value(fine, mask);
                        return fine;
                    }
        };

    template <> struct Prolongation<tags::GPU::CUDA, NONE>
    {
        static DenseVector<float> & value(DenseVector<float> & fine,
                const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask, HONEI_UNUSED BandedMatrixQ1<float> & prolmat);

        static DenseVector<double> & value(DenseVector<double> & fine,
                const DenseVector<double> & coarse, const DenseVector<unsigned long> & mask, HONEI_UNUSED BandedMatrixQ1<double> & prolmat);

        static DenseVector<float> & value(DenseVector<float> & fine,
                const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask, HONEI_UNUSED SparseMatrixELL<float> & prolmat);

        static DenseVector<double> & value(DenseVector<double> & fine,
                const DenseVector<double> & coarse, const DenseVector<unsigned long> & mask, HONEI_UNUSED SparseMatrixELL<double> & prolmat);
    };
}
#endif
