/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. Math is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * Math is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#pragma once
#ifndef LIBMATH_GUARD_DEFECT_HH
#define LIBMATH_GUARD_DEFECT_HH 1

#include<honei/la/banded_matrix_qx.hh>
#include<honei/la/sparse_matrix_ell.hh>
#include<honei/la/dense_vector.hh>
#include<honei/la/algorithm.hh>
#include<honei/la/product.hh>
#include<honei/la/difference.hh>
#include<honei/backends/multicore/thread_pool.hh>
#include<honei/util/operation_wrapper.hh>
#include<honei/util/profiler.hh>
using namespace honei;
namespace honei
{
    template<typename Tag_ = tags::CPU> struct Defect;

    template <> struct Defect<tags::CPU>
    {
        public:
            template<typename DT_>
                static DenseVector<DT_> value(const DenseVector<DT_> & right_hand_side, const BandedMatrixQx<Q1Type, DT_> & system, const DenseVector<DT_> & x)
                {
                    DenseVector<DT_> result(right_hand_side.copy());
                    DenseVector<DT_> prod(Product<tags::CPU>::value(system, x));
                    Difference<tags::CPU>::value(result,prod);
                    return result;

                    /*if (x.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(x.size(), system.columns());
                    }
                    if (right_hand_side.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                    }

                    right_hand_side.lock(lm_read_only);
                    system.lock(lm_read_only);
                    x.lock(lm_read_only);
                    DenseVector<DT_> result(right_hand_side.size());
                    result.lock(lm_write_only);
                    unsigned long n = right_hand_side.size();
                    unsigned long root_n = (unsigned long)sqrt(n);

                    DT_ * rhs = right_hand_side.elements();
                    DT_ * x_old = x.elements();
                    DT_ * x_new = result.elements();

                    DT_ * ll = system.band(LL).elements();
                    DT_ * ld = system.band(LD).elements();
                    DT_ * lu = system.band(LU).elements();

                    DT_ * dl = system.band(DL).elements();
                    DT_ * dd = system.band(DD).elements();
                    DT_ * du = system.band(DU).elements();

                    DT_ * ul = system.band(UL).elements();
                    DT_ * ud = system.band(UD).elements();
                    DT_ * uu = system.band(UU).elements();

                    unsigned long i(0);
                    //index 0
                    x_new[i] = ((rhs[i] - (dd[i] * x_old[i] +
                                    du[i] * x_old[1] +
                                    ul[i] * x_old[root_n - 1] +
                                    ud[i] * x_old[root_n] +
                                    uu[i] * x_old[root_n + 1])));

                    //index in [1, root_n -1[
                    i = 1;
                    for(; i < root_n - 1 ; ++i)
                    {
                        x_new[i] = ((rhs[i] - (dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[i + 1] +
                                        ul[i] * x_old[i + root_n - 1] +
                                        ud[i] * x_old[i + root_n] +
                                        uu[i] * x_old[i + root_n + 1])));
                    }

                    //index root_n -1
                    i = root_n - 1;
                    x_new[i] = ((rhs[i] - (lu[i] * x_old[i - (root_n - 1)] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n] +
                                    uu[i] * x_old[i + root_n + 1])));

                    //index root_n
                    i = root_n;
                    x_new[i] = ((rhs[i] - (ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - (root_n - 1)] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n] +
                                    uu[i] * x_old[i + root_n + 1])));

                    //index in [root_n + 1, n - (root_n + 1)[
                    i = root_n + 1;
                    for(; i < n - (root_n  + 1) ; ++i)
                    {
                        x_new[i] = ((rhs[i] - (ll[i] * x_old[i - root_n - 1] +
                                        ld[i] * x_old[i - root_n] +
                                        lu[i] * x_old[i - root_n + 1] +
                                        dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[ i + 1] +
                                        ul[i] * x_old[i + root_n - 1] +
                                        ud[i] * x_old[i + root_n] +
                                        uu[i] * x_old[i + root_n + 1])));
                    }

                    //index n - (root_n + 1)
                    i = n - (root_n + 1);
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n])));

                    //index n - root_n
                    i = n - root_n;
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1])));

                    //index in [n - root_n + 1, n -1[
                    i = n - root_n + 1;
                    for(; i < n - 1; ++i)
                    {
                        x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                        ld[i] * x_old[i - root_n] +
                                        lu[i] * x_old[i - root_n + 1] +
                                        dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[ i + 1])));
                    }

                    //index n - 1
                    i = n - 1;
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i])));

                    right_hand_side.unlock(lm_read_only);
                    system.unlock(lm_read_only);
                    x.unlock(lm_read_only);
                    result.unlock(lm_write_only);
                    return result;
                    */
                }

            template<typename DT_>
                static DenseVector<DT_> & value(DenseVector<DT_> & result, const DenseVector<DT_> & right_hand_side, const BandedMatrixQx<Q1Type, DT_> & system, const DenseVector<DT_> & x)
                {
                    //DenseVector<DT_> result(right_hand_side.copy());
                    Difference<tags::CPU>::value(result, right_hand_side, Product<tags::CPU>::value(system, x) );
                    return result;
                    /*
                    if (x.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(x.size(), system.columns());
                    }
                    if (right_hand_side.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                    }

                    right_hand_side.lock(lm_read_only);
                    system.lock(lm_read_only);
                    x.lock(lm_read_only);
                    result.lock(lm_write_only);
                    unsigned long n = right_hand_side.size();
                    unsigned long root_n = (unsigned long)sqrt(n);

                    DT_ * rhs = right_hand_side.elements();
                    DT_ * x_old = x.elements();
                    DT_ * x_new = result.elements();

                    DT_ * ll = system.band(LL).elements();
                    DT_ * ld = system.band(LD).elements();
                    DT_ * lu = system.band(LU).elements();

                    DT_ * dl = system.band(DL).elements();
                    DT_ * dd = system.band(DD).elements();
                    DT_ * du = system.band(DU).elements();

                    DT_ * ul = system.band(UL).elements();
                    DT_ * ud = system.band(UD).elements();
                    DT_ * uu = system.band(UU).elements();

                    unsigned long i(0);
                    //index 0
                    x_new[i] = ((rhs[i] - (dd[i] * x_old[i] +
                                    du[i] * x_old[1] +
                                    ul[i] * x_old[root_n - 1] +
                                    ud[i] * x_old[root_n] +
                                    uu[i] * x_old[root_n + 1])));

                    //index in [1, root_n -1[
                    i = 1;
                    for(; i < root_n - 1 ; ++i)
                    {
                        x_new[i] = ((rhs[i] - (dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[i + 1] +
                                        ul[i] * x_old[i + root_n - 1] +
                                        ud[i] * x_old[i + root_n] +
                                        uu[i] * x_old[i + root_n + 1])));
                    }

                    //index root_n -1
                    i = root_n - 1;
                    x_new[i] = ((rhs[i] - (lu[i] * x_old[i - (root_n - 1)] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n] +
                                    uu[i] * x_old[i + root_n + 1])));

                    //index root_n
                    i = root_n;
                    x_new[i] = ((rhs[i] - (ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - (root_n - 1)] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n] +
                                    uu[i] * x_old[i + root_n + 1])));

                    //index in [root_n + 1, n - (root_n + 1)[
                    i = root_n + 1;
                    for(; i < n - (root_n  + 1) ; ++i)
                    {
                        x_new[i] = ((rhs[i] - (ll[i] * x_old[i - root_n - 1] +
                                        ld[i] * x_old[i - root_n] +
                                        lu[i] * x_old[i - root_n + 1] +
                                        dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[ i + 1] +
                                        ul[i] * x_old[i + root_n - 1] +
                                        ud[i] * x_old[i + root_n] +
                                        uu[i] * x_old[i + root_n + 1])));
                    }

                    //index n - (root_n + 1)
                    i = n - (root_n + 1);
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1] +
                                    ud[i] * x_old[i + root_n])));

                    //index n - root_n
                    i = n - root_n;
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i] +
                                    du[i] * x_old[ i + 1] +
                                    ul[i] * x_old[i + root_n - 1])));

                    //index in [n - root_n + 1, n -1[
                    i = n - root_n + 1;
                    for(; i < n - 1; ++i)
                    {
                        x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                        ld[i] * x_old[i - root_n] +
                                        lu[i] * x_old[i - root_n + 1] +
                                        dl[i] * x_old[i - 1] +
                                        dd[i] * x_old[i] +
                                        du[i] * x_old[ i + 1])));
                    }

                    //index n - 1
                    i = n - 1;
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                                    ld[i] * x_old[i - root_n] +
                                    lu[i] * x_old[i - root_n + 1] +
                                    dl[i] * x_old[i - 1] +
                                    dd[i] * x_old[i])));

                    right_hand_side.unlock(lm_read_only);
                    system.unlock(lm_read_only);
                    x.unlock(lm_read_only);
                    result.unlock(lm_write_only);
                    return result;
                    */
                }

            template<typename DT_>
                static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                {
                    if (x.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(x.size(), system.columns());
                    }
                    if (right_hand_side.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                    }

                    DenseVector<DT_> result(right_hand_side.copy());
                    DenseVector<DT_> temp(right_hand_side.size());
                    Product<tags::CPU>::value(temp, system, x);
                    Difference<tags::CPU>::value(result, temp);
                    return result;
                }

            template<typename DT_>
                static DenseVector<DT_> & value(DenseVector<DT_> & result, DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                {
                    if (x.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(x.size(), system.columns());
                    }
                    if (right_hand_side.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                    }

                    DenseVector<DT_> temp(right_hand_side.size());
                    Product<tags::CPU>::value(temp, system, x);

                    result.lock(lm_write_only);
                    right_hand_side.lock(lm_read_only);
                    temp.lock(lm_read_only);
                    for(unsigned long i(0) ; i < result.size() ; ++i)
                    {
                        result[i] = right_hand_side[i] - temp[i];
                    }
                    result.unlock(lm_write_only);
                    right_hand_side.unlock(lm_read_only);
                    temp.unlock(lm_read_only);

                    return result;
                }

            /*template <typename DT_>
              static inline BenchmarkInfo get_benchmark_info(const DenseVector<DT_> & right_hand_side, const SparseMatrixELL<DT_> & system, const DenseVector<DT_> & x)
              {
              BenchmarkInfo result;
              result.flops = :;
              result.load = a.rows() * a.columns() * (sizeof(DT1_) + sizeof(DT2_));
              result.store = a.rows() * a.columns() * sizeof(DT1_);
              result.size.push_back(a.rows() * a.columns());
              result.size.push_back(b.rows() * b.columns());
              return result;
              }*/
    };

    template<>
        struct Defect<tags::GPU::CUDA>
        {
            public:
                static DenseVector<float> value(const DenseVectorContinuousBase<float> & right_hand_side,
                        const BandedMatrixQx<Q1Type, float> & system, const DenseVectorContinuousBase<float> & x);

                static DenseVector<double> value(const DenseVectorContinuousBase<double> & right_hand_side,
                        const BandedMatrixQx<Q1Type, double> & system, const DenseVectorContinuousBase<double> & x);

                static DenseVector<float> & value(DenseVector<float> & result, const DenseVectorContinuousBase<float> & right_hand_side,
                        const BandedMatrixQx<Q1Type, float> & system, const DenseVectorContinuousBase<float> & x);

                static DenseVector<double> & value(DenseVector<double> & result, const DenseVectorContinuousBase<double> & right_hand_side,
                        const BandedMatrixQx<Q1Type, double> & system, const DenseVectorContinuousBase<double> & x);

                static DenseVector<float> value(const DenseVectorContinuousBase<float> & right_hand_side,
                        const SparseMatrixELL<float> & system, const DenseVectorContinuousBase<float> & x);

                static DenseVector<double> value(const DenseVectorContinuousBase<double> & right_hand_side,
                        const SparseMatrixELL<double> & system, const DenseVectorContinuousBase<double> & x);

                static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & right_hand_side,
                        const SparseMatrixELL<float> & system, const DenseVectorContinuousBase<float> & x);

                static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & right_hand_side,
                        const SparseMatrixELL<double> & system, const DenseVectorContinuousBase<double> & x);

        };

    template<>
        struct Defect<tags::OpenCL::CPU>
        {
            public:
                template <typename DT_>
                static DenseVector<DT_> value(const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const BandedMatrixQx<Q1Type, DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const BandedMatrixQx<Q1Type, DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVector<DT_> value(const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const SparseMatrixELL<DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const SparseMatrixELL<DT_> & system, const DenseVectorContinuousBase<DT_> & x);
        };

    template<>
        struct Defect<tags::OpenCL::GPU>
        {
            public:
                template <typename DT_>
                static DenseVector<DT_> value(const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const BandedMatrixQx<Q1Type, DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const BandedMatrixQx<Q1Type, DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVector<DT_> value(const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const SparseMatrixELL<DT_> & system, const DenseVectorContinuousBase<DT_> & x);

                template <typename DT_>
                static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const DenseVectorContinuousBase<DT_> & right_hand_side,
                        const SparseMatrixELL<DT_> & system, const DenseVectorContinuousBase<DT_> & x);
        };

    template<>
        struct Defect<tags::CPU::SSE>
        {

            public:
                static DenseVector<double> value(const DenseVector<double> & right_hand_side, const BandedMatrixQx<Q1Type, double> & system, const DenseVector<double> & x);
                static DenseVector<float> value(const DenseVector<float> & right_hand_side, const BandedMatrixQx<Q1Type, float> & system, const DenseVector<float> & x);

                static DenseVector<double> & value(DenseVector<double> & result, const DenseVector<double> & right_hand_side, const BandedMatrixQx<Q1Type, double> & system, const DenseVector<double> & x);
                static DenseVector<float> & value(DenseVector<float> & result, const DenseVector<float> & right_hand_side, const BandedMatrixQx<Q1Type, float> & system, const DenseVector<float> & x);

                template<typename DT_>
                    static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                    {
                        PROFILER_START("Defect SMELL tags::CPU::SSE");
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> result(right_hand_side.size());
                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::CPU::SSE>::value(temp, system, x);
                        Difference<tags::CPU::SSE>::value(result, right_hand_side, temp);
                        PROFILER_STOP("Defect SMELL tags::CPU::SSE");
                        return result;
                    }

                    static DenseVector<float> & value(DenseVector<float> & result, const DenseVector<float> & right_hand_side, const SparseMatrixELL<float> & system, const DenseVector<float> & x, unsigned long row_start = 0, unsigned long row_end = 0);

                    static DenseVector<double> & value(DenseVector<double> & result, const DenseVector<double> & right_hand_side, const SparseMatrixELL<double> & system, const DenseVector<double> & x, unsigned long row_start = 0, unsigned long row_end = 0);
        };

    template<>
        struct Defect<tags::GPU::MultiCore::CUDA>
        {

            //todo DIRK
            //TODO use atomar MC cuda defect
            public:
                template<typename DT_>
                    static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> result(right_hand_side.size());
                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::GPU::MultiCore::CUDA>::value(temp, system, x);
                        Difference<tags::GPU::MultiCore::CUDA>::value(result, right_hand_side, temp);
                        return result;
                    }

                template<typename DT_>
                    static DenseVector<DT_> & value(DenseVector<DT_> & result, DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::GPU::MultiCore::CUDA>::value(temp, system, x);
                        Difference<tags::GPU::MultiCore::CUDA>::value(result, right_hand_side, temp);

                        return result;
                    }
        };

    template<>
        struct Defect<tags::CPU::MultiCore>
        {

            public:
                template<typename DT_>
                    static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, BandedMatrixQx<Q1Type, DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> result(right_hand_side.size());
                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::CPU::MultiCore>::value(temp, system, x);
                        Difference<tags::CPU::MultiCore>::value(result, right_hand_side, temp);
                        return result;
                    }

                template<typename DT_>
                    static DenseVector<DT_> & value(DenseVector<DT_> & result, DenseVector<DT_> & right_hand_side, BandedMatrixQx<Q1Type, DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::CPU::MultiCore>::value(temp, system, x);
                        Difference<tags::CPU::MultiCore>::value(result, right_hand_side, temp);

                        return result;
                    }

                template<typename DT_>
                    static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> result(right_hand_side.size());
                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::CPU::MultiCore>::value(temp, system, x);
                        Difference<tags::CPU::MultiCore>::value(result, right_hand_side, temp);
                        return result;
                    }

                template<typename DT_>
                    static DenseVector<DT_> & value(DenseVector<DT_> & result, DenseVector<DT_> & right_hand_side, SparseMatrixELL<DT_> & system, DenseVector<DT_> & x)
                    {
                        if (x.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(x.size(), system.columns());
                        }
                        if (right_hand_side.size() != system.columns())
                        {
                            throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                        }

                        DenseVector<DT_> temp(right_hand_side.size());
                        Product<tags::CPU::MultiCore>::value(temp, system, x);
                        Difference<tags::CPU::MultiCore>::value(result, right_hand_side, temp);

                        return result;
                    }
        };

    namespace mc
    {
        template <typename Tag_> struct Defect
        {

                template<typename DT_>
                static DenseVector<DT_> & value(DenseVector<DT_> & result, DenseVector<DT_> & right_hand_side, BandedMatrixQx<Q1Type, DT_> & system, DenseVector<DT_> & x)
                {
                    if (x.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(x.size(), system.columns());
                    }
                    if (right_hand_side.size() != system.columns())
                    {
                        throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                    }

                    DenseVector<DT_> temp(right_hand_side.size());
                    Product<Tag_>::value(temp, system, x);
                    Difference<Tag_>::value(result, right_hand_side, temp);

                    return result;
                }

                template <typename DT_>
                    static DenseVector<DT_> value(const DenseVector<DT_> & rhs, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
                    {
                        if (b.size() != a.columns())
                        {
                            throw VectorSizeDoesNotMatch(b.size(), a.columns());
                        }
                        if (a.rows() != rhs.size())
                        {
                            throw VectorSizeDoesNotMatch(a.rows(), rhs.size());
                        }

                        DenseVector<DT_> result(a.rows());
                        //fill<typename Tag_::DelegateTo>(result, DT_(0));

                        unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                                    mc::ThreadPool::instance()->num_threads()));

                        TicketVector tickets;

                        unsigned long limits[max_count + 1];
                        limits[0] = 0;
                        for (unsigned long i(1) ; i < max_count; ++i)
                        {
                            limits[i] = limits[i-1] + a.rows() / max_count;
                        }
                        limits[max_count] = a.rows();

                        for (unsigned long i(0) ; i < max_count ; ++i)
                        {
                            OperationWrapper<honei::Defect<typename Tag_::DelegateTo>, DenseVector<DT_>, DenseVector<DT_>,
                                DenseVector<DT_>, SparseMatrixELL<DT_>, DenseVector<DT_>, unsigned long, unsigned long > wrapper(result);
                            tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, result, rhs, a, b, limits[i], limits[i+1])));
                        }

                        tickets.wait();

                        return result;
                    }

                template <typename DT_>
                    static DenseVector<DT_> value(DenseVector<DT_> & result, const DenseVector<DT_> & rhs, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
                    {
                        if (b.size() != a.columns())
                        {
                            throw VectorSizeDoesNotMatch(b.size(), a.columns());
                        }
                        if (a.rows() != result.size())
                        {
                            throw VectorSizeDoesNotMatch(a.rows(), result.size());
                        }

                        //fill<typename Tag_::DelegateTo>(result, DT_(0));

                        unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                                    mc::ThreadPool::instance()->num_threads()));

                        TicketVector tickets;

                        unsigned long limits[max_count + 1];
                        limits[0] = 0;
                        for (unsigned long i(1) ; i < max_count; ++i)
                        {
                            limits[i] = limits[i-1] + a.rows() / max_count;
                        }
                        limits[max_count] = a.rows();

                        for (unsigned long i(0) ; i < max_count ; ++i)
                        {
                            OperationWrapper<honei::Defect<typename Tag_::DelegateTo>, DenseVector<DT_>, DenseVector<DT_>,
                                DenseVector<DT_>, SparseMatrixELL<DT_>, DenseVector<DT_>, unsigned long, unsigned long > wrapper(result);
                            tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, result, rhs, a, b, limits[i], limits[i+1])));
                        }

                        tickets.wait();

                        return result;
                    }
        };

    }

    template <> struct Defect<tags::CPU::MultiCore::SSE> :
        public mc::Defect<tags::CPU::MultiCore::SSE>
        {
        };
}

#endif
