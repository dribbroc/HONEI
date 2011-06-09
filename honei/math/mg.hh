/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
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


#ifndef MATH_GUARD_MG_HH
#define MATH_GUARD_MG_HH 1

#include <operator.hh>

namespace honei
{
    template<typename Tag_, typename NormType_>
    struct MG
    {
        public:
            template<typename MatrixType_, typename VectorType_, typename DT_>
            static void value(MatrixType_ & A, VectorType_ & b, VectorType_ & x, OperatorList & cycle, unsigned long max_iters, unsigned long & used_iters, DT_ eps_relative)
            {
                VectorType_ r(x.size());
                Defect<Tag_>::value(r, b, A, x);
                DT_ rnorm_initial(NormType_::value(r));
                DT_ rnorm_current(1e16);

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                    Defect<Tag_>::value(r, b, A, x);
                    rnorm_current = NormType_::value(r);

                    used_iters = i + 1;

                    if(rnorm_current < eps_relative * rnorm_initial)
                        break;
                }
            }
    };

    struct MGSmoother
    {
        public:
            static void value(OperatorList & cycle, unsigned long max_iters)
            {
                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                }
            }
    };
}

#endif
