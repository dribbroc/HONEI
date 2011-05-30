/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
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

#pragma once
#ifndef LIBMATH_GUARD_RICHARDSON_HH
#define LIBMATH_GUARD_RICHARDSON_HH 1

#include <honei/util/tags.hh>
#include <honei/la/product.hh>
#include <honei/la/sum.hh>
#include <honei/la/scale.hh>
#include <honei/la/norm.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/la/algorithm.hh>
#include <honei/math/defect.hh>
#include <honei/math/preconditioning.hh>

namespace honei
{
    template <typename MatrixSourceType_, typename MatrixTargetType_, typename DataType_>
    class LSData
    {
        private:
            LSData()
            {
            }

        public:
            LSData(MatrixSourceType_ * A,
                   MatrixTargetType_ * C,
                   DenseVector<DataType_> * b,
                   DenseVector<DataType_> * y,
                   DenseVector<DataType_> * t0,
                   DenseVector<DataType_> * t1) :
                system_matrix(A),
                precon_matrix(C),
                rhs(b),
                result(y),
                temp_0(t0),
                temp_1(t1)
            {
            }
            MatrixSourceType_ * system_matrix;
            MatrixTargetType_ * precon_matrix;
            DenseVector<DataType_> * rhs;
            DenseVector<DataType_> * result;
            DenseVector<DataType_> * temp_0;
            DenseVector<DataType_> * temp_1;
    };

    class LSInfo
    {
        private:
            /*LSInfo() :
            {
            }*/

        public:
            LSInfo(bool ccr,
                   bool d,
                   double tr,
                   double df,
                   unsigned long mi,
                   unsigned long & i) :
                convergence_control_relative(ccr),
                damping(d),
                tolerance_relative(tr),
                damping_factor(df),
                max_iters(mi),
                iters(i)
            {
            }
            bool convergence_control_relative;
            bool damping;
            double tolerance_relative;
            double damping_factor;
            unsigned long max_iters;
            unsigned long & iters;
    };

    /**
     * \brief Solution of linear system with Richardson.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    template <typename Tag_ = tags::CPU, typename Preconditioner_ = Preconditioning<Tag_, methods::NONE> >
    struct Richardson;

    template <typename Tag_, typename Preconditioner_>
    struct Richardson
    {
        public:
            /**
            * \brief Returns solution of linear system with the Richardson method.
            *
            * \param data The linear system.
            * \param config Solver configuration.
            *
            */
            template <typename MatrixSourceType_, typename MatrixTargetType_, typename DataType_>
            static void value(
                    LSData<MatrixSourceType_, MatrixTargetType_, DataType_> & data,
                    LSInfo & info)
            {
                CONTEXT("When solving linear system with Richardson: ");

                Preconditioner_::value(*data.precon_matrix, *data.system_matrix);
                if(info.damping)
                {
                    Scale<Tag_>::value(data.precon_matrix->Ax(), DataType_(info.damping_factor));
                }

                DataType_ initial_defect_norm(1000);
                if(info.convergence_control_relative)
                {
                    Defect<Tag_>::value(*data.temp_0, *data.rhs, *data.system_matrix, *data.result);
                    initial_defect_norm = Norm<vnt_l_two, false, Tag_>::value(*data.temp_0);
                }
                for(unsigned long i(0) ; i < info.max_iters ; ++i)
                {
                    Defect<Tag_>::value(*data.temp_0, *data.rhs, *data.system_matrix, *data.result);
                    Product<Tag_>::value(*data.temp_1, *data.precon_matrix, *data.temp_0);
                    Sum<Tag_>::value(*data.result, *data.temp_1);

                    if(info.convergence_control_relative)
                    {
                        //Defect<Tag_>::value(*data.temp_0, *data.rhs, *data.system_matrix, *data.result);
                        if(Norm<vnt_l_two, false, Tag_>::value(*data.temp_0) < info.tolerance_relative * initial_defect_norm)
                        {
                            info.iters = i;
                            break;
                        }
                    }
                    if(i == info.max_iters - 1)
                        info.iters = i;
                }
            }
    };
}
#endif
