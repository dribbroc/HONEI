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

#pragma once
#ifndef MATH_GUARD_IR_HH
#define MATH_GUARD_IR_HH 1

#include<honei/math/mg.hh>

//TODO: plug in solver as type
//TODO: innertag is only implicitly given in OperatorList cycle
namespace honei
{
    template<typename OuterTag_, typename InnerTag_, typename NormType_>
    struct IRSolver
    {
        public:
            template<typename OuterMatrixType_,
                     typename OuterVectorType_,
                     typename InnerMatrixType_,
                     typename InnerVectorType_,
                     typename InnerTransferContType_,
                     typename InnerPreconContType_>
            static void value(OuterMatrixType_ A, OuterVectorType_ b, OuterVectorType_ x, MGData<InnerMatrixType_, InnerVectorType_, InnerTransferContType_, InnerPreconContType_> & data, OperatorList & cycle, double eps, unsigned long max_iters, unsigned long & used_iters)
            {
                CONTEXT("When solving linear system with Iterative Refinement :");
                ASSERT(cycle.size() > 0, "OperatorList is empty!");
                PROFILER_START("IterativeRefinementMGSolver");

                OuterVectorType_ r(x.size());//TODO:move allocation to preprocessing
                OuterVectorType_ c(x.size());//TODO:move allocation to preprocessing
                double rnorm_initial(1e16);
                double rnorm_current(1e16);

                //std::cout << "starting cycles" << std::endl;
                //inner MG-cycles, TODO: see above
                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    Defect<OuterTag_>::value(r, b, A, x);
                    if(i == 0)
                        rnorm_initial = NormType_::value(r);
                    else
                    {
                        rnorm_current = NormType_::value(r);
                        used_iters = i;

                        if(rnorm_current < eps * rnorm_initial)
                            break;
                    }

                    convert<OuterTag_>(data.b.at(data.b.size() - 1), r);//TODO: target, then source?

                    MGSolver<InnerTag_, NormType_>::value(data, cycle);

                    convert<OuterTag_>(c, data.x.at(data.x.size() - 1));//TODO: target, then source?

                    Sum<OuterTag_>::value(x, c);

                }

                PROFILER_STOP("IterativeRefinementMGSolver");
            }
    };
}
#endif
