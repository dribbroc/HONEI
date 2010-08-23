/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2009 Sven Mallach <mallach@honei.org>
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
             typename T4_ = void, typename T5_ = void, typename T6_ = void>
    struct OperationWrapper;

    template <typename Operation_, typename ResultType_, typename T1_>
    struct OperationWrapper<Operation_, ResultType_, T1_>
    {
        typedef void result_type;
        ResultType_ & result;

        OperationWrapper(ResultType_ & r) :
            result(r)
        {
        }

        void operator() (T1_ & param_one)
        {
            result = Operation_::value(param_one);
        }
    };

    template <typename Operation_, typename ResultType_, typename T1_, typename T2_>
    struct OperationWrapper<Operation_, ResultType_, T1_, T2_>
    {
        typedef void result_type;
        ResultType_ & result;

        OperationWrapper(ResultType_ & r) :
            result(r)
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
        typedef void result_type;
        ResultType_ & result;

        OperationWrapper(ResultType_ & r) :
            result(r)
        {
        }

        void operator() (T1_ & param_one, const T2_ & param_two, const T3_ & param_three)
        {
            result = Operation_::value(param_one, param_two, param_three);
        }
    };

    template <typename Operation_, typename ResultType_, typename T1_, typename T2_, typename T3_, typename T4_>
    struct OperationWrapper<Operation_, ResultType_, T1_, T2_, T3_, T4_>
    {
        typedef void result_type;
        ResultType_ & result;

        OperationWrapper(ResultType_ & r) :
            result(r)
        {
        }

        void operator() (T1_ & param_one, const T2_ & param_two, const T3_ & param_three, const T4_ & param_four)
        {
            result = Operation_::value(param_one, param_two, param_three, param_four);
        }
    };

    template <typename Operation_, typename ResultType_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_>
    struct OperationWrapper<Operation_, ResultType_, T1_, T2_, T3_, T4_, T5_>
    {
        typedef void result_type;
        ResultType_ & result;

        OperationWrapper(ResultType_ & r) :
            result(r)
        {
        }

        void operator() (T1_ & param_one, const T2_ & param_two, const T3_ & param_three, const T4_ & param_four, const T5_ & param_five)
        {
            result = Operation_::value(param_one, param_two, param_three, param_four, param_five);
        }
    };
}

#endif
