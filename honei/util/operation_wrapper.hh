/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef HONEI_GUARD_UTIL_OPERATION_WRAPPER_HH
#define HONEI_GUARD_UTIL_OPERATION_WRAPPER_HH 1

namespace honei
{
    /**
     * OperationWrapper is used to construct a bindable functor for calls to Operations.
     *
     * The reasoning for a custom wrapper is that virtually all operations have
     * overloaded functions.
     */
    template <typename Operation_, typename ResultType_, typename T1_, typename T2_ = void, typename T3_ = void,
             typename T4_ = void, typename T4_ = void, typename T5_ = void, typename T6_ = void>
    struct OperationWrapper;

    template <typename Operation_, typename ResultType_, typename T1_>
    struct OperationWrapper<Operation_, ResultType_, T1_>
    {
        ResultType_ result;

        OperationWrapper(ResultType_ & result) :
            result(result)
        {
        }

        void operator() (T1_ & param_one)
        {
            result = Operation_::value(parame_one);
        }
    };

    template <typename Operation_, typename ResultType_, typename T1_, typename T2_>
    struct OperationWrapper<Operation_, ResultType_, T1_, T2_>
    {
        ResultType_ result;

        OperationWrapper(ResultType_ & result) :
            result(result)
        {
        }

        void operator() (T1_ & param_one, const T2_ & param_two)
        {
            result = Operation_::value(param_one, param_two);
        }
    };

    template <typename Operation_, typename ResultType_, typename T1_, typename T2_, typename T3_>
    struct OperationWrapper<Operation_, ResultType_, T1_, T2_, T3_>
    {
        ResultType_ result;

        OperationWrapper(ResultType_ & result) :
            result(result)
        {
        }

        void operator() (T1_ & param_one, const T2_ & param_two, const T3_ & param_three)
        {
            result = Operation_::value(param_one, param_two, param_three);
        }
    };
}

#endif
