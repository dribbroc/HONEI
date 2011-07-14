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
#ifndef MATH_GUARD_OPERATOR_HH
#define MATH_GUARD_OPERATOR_HH 1

#include <honei/math/defect.hh>
#include <honei/la/norm.hh>
#include <honei/la/sum.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/math/cg.hh>
#include <honei/math/ri.hh>

namespace honei
{
    class Operator
    {
        public:
            virtual void value() = 0;

            virtual ~Operator()
            {
            }

            virtual std::string to_string() = 0;
    };

    template<typename Tag_, typename VectorType_>
    class CopyOperator : public Operator
    {
        public:

            virtual std::string to_string()
            {
                std::string result("CopyOperator\n");
                return result;;
            }

            CopyOperator(VectorType_& x, VectorType_& y) :
                _x(x),
                _y(y)
        {
            CONTEXT("When creating CopyOperator:");
        }

            virtual void value()
            {
                CONTEXT("When evaluating CopyOperator:");
                copy<Tag_>(_x, _y);
            }

            virtual ~CopyOperator()
            {
                CONTEXT("When destroying CopyOperator:");
            }

        private:
            VectorType_ _x;
            VectorType_ _y;
    };

    template<typename Tag_, typename MatType_, typename VectorType_>
    class DefectOperator : public Operator
    {
        public:
            virtual std::string to_string()
            {
                return "DefectOperator";
            }

            DefectOperator(VectorType_& y, VectorType_& b, MatType_& A, VectorType_& x) :
                _y(y),
                _b(b),
                _A(A),
                _x(x)
            {
                CONTEXT("When creating DefectOperator:");
            }

            virtual void value()
            {
                CONTEXT("When evaluating DefectOperator:");
                Defect<Tag_>::value(_y, _b, _A, _x);
                //std::cout << _y << std::endl;
            }

            virtual ~DefectOperator()
            {
                CONTEXT("When destroying DefectOperator:");
            }

        private:
            VectorType_ _y;
            VectorType_ _b;
            MatType_ _A;
            VectorType_ _x;
    };

    template<typename Tag_, VectorNormType normtype, typename DT_, bool root>
    class NormOperator : public Operator
    {
        public:
            virtual std::string to_string()
            {
                return "NormOperator";
            }

            NormOperator(DT_& y, DenseVector<DT_>& x) :
                _x(x),
                _y(y)
        {
            CONTEXT("When creating NormOperator:");
        }

            virtual void value()
            {
                CONTEXT("When evaluating NormOperator:");
                _y = Norm<normtype, root, Tag_>::value(_x);
            }

            virtual ~NormOperator()
            {
                CONTEXT("When destroying NormOperator:");
            }

        private:
            DenseVector<DT_> _x;
            DT_& _y;
    };

    template<typename Tag_, typename VectorType_>
    class SumOperator : public Operator
    {
        public:

            virtual std::string to_string()
            {
                std::string result("SumOperator\n");
                result += "  left.size()= ";
                result += stringify(_x.size());
                result += "\n";
                result += "  right.size()= ";
                result += stringify(_y.size());
                result += "\n";
                return result;;
            }

            SumOperator(VectorType_& y, VectorType_& x) :
                _x(x),
                _y(y)
        {
            CONTEXT("When creating SumOperator:");
        }

            virtual void value()
            {
                CONTEXT("When evaluating SumOperator:");
                Sum<Tag_>::value(_y, _x);
            }

            virtual ~SumOperator()
            {
                CONTEXT("When destroying SumOperator:");
            }

        private:
            VectorType_ _x;
            VectorType_ _y;
    };

    template<typename Tag_, typename DT_>
    class ScaledSumOperator : public Operator
    {
        public:
            virtual std::string to_string()
            {
                return "ScaledSumOperator";
            }

            ScaledSumOperator(DenseVector<DT_>& y, DenseVector<DT_>& x, DT_ alpha) :
                _x(x),
                _y(y),
                _alpha(alpha)
        {
            CONTEXT("When creating ScaledSumOperator:");
        }

            virtual void value()
            {
                CONTEXT("When evaluating ScaledSumOperator:");
                ScaledSum<Tag_>::value(_y, _x, _alpha);
            }

            virtual ~ScaledSumOperator()
            {
                CONTEXT("When destroying ScaledSumOperator:");
            }

        private:
            DenseVector<DT_> _x;
            DenseVector<DT_> _y;
            DT_ _alpha;
    };

    template<typename SolverType_, typename MatrixType_, typename VectorType_>
    class SolverOperator : public Operator
    {
        //TODO: how to build in preconditioning -> standard-value?
        public:
            virtual std::string to_string()
            {
                return "SolverOperator";
            }

            SolverOperator(MatrixType_ & A,
                           VectorType_ & b,
                           VectorType_ & x,
                           unsigned long max_iters,
                           unsigned long & used_iters,
                           double eps_relative = 1e-8) :
                _A(A),
                _b(b),
                _x(x),
                _max_iters(max_iters),
                _used_iters(used_iters),
                _eps_relative(eps_relative)
            {
                CONTEXT("When creating SolverOperator:");
            }

            virtual void value()
            {
                CONTEXT("When evaluating SolverOperator:");
                //std::cout << _x << std::endl;
                SolverType_::value(_A, _b, _x, _max_iters, _used_iters, _eps_relative);
                //std::cout << _x << std::endl;
            }

            virtual ~SolverOperator()
            {
                CONTEXT("When destroying SolverOperator:");
            }

        private:
            MatrixType_ _A;
            VectorType_ _b;
            VectorType_ _x;
            unsigned long _max_iters;
            unsigned long & _used_iters;
            double _eps_relative;
    };

    template<typename SmootherType_, typename MatrixType_, typename VectorType_, typename PreconContType_>
    class SmootherOperator : public Operator
    {
        //TODO: more temps, non-matrix preconditioning
        public:
            virtual std::string to_string()
            {
                return "SmootherOperator";
            }

            SmootherOperator(MatrixType_ & A,
                             PreconContType_ & P,
                             VectorType_ & b,
                             VectorType_ & x,
                             VectorType_ & temp_0,
                             VectorType_ & temp_1,
                             unsigned long max_iters) :
                _A(A),
                _P(P),
                _b(b),
                _x(x),
                _temp_0(temp_0),
                _temp_1(temp_1),
                _max_iters(max_iters)
            {
                CONTEXT("When creating SmootherOperator:");
            }

            virtual void value()
            {
                CONTEXT("When evaluating SmootherOperator:");
                std::cout << "x before" << _x;
                std::cout << "rhs" << _b;
                SmootherType_::value(_A, _P, _b, _x, _temp_0, _temp_1, _max_iters);
                std::cout << "x after" << _x;
            }

            virtual ~SmootherOperator()
            {
                CONTEXT("When destroying SmootherOperator:");
            }

        private:
            MatrixType_ _A;
            PreconContType_ _P;
            VectorType_ _b;
            VectorType_ _x;
            VectorType_ _temp_0;
            VectorType_ _temp_1;
            unsigned long _max_iters;
    };

    template<typename TransferType_, typename MatrixType_, typename VectorType_>
    class TransferOperator : public Operator
    {
        public:
            virtual std::string to_string()
            {
                std::string result("TransferOperator\n");
                result += "  left.size()= ";
                result += stringify(_left.size());
                result += "\n";
                result += "  right.size()= ";
                result += stringify(_right.size());
                result += "\n";

                result += "  mat_dims= ";
                result += stringify(_mat.rows());
                result += " x ";
                result += stringify(_mat.columns());
                result += "\n";

                return result;
            }

            TransferOperator(VectorType_ & left, VectorType_ & right, MatrixType_ & mat) :
                _left(left),
                _right(right),
                _mat(mat)
            {
                CONTEXT("When creating TransferOperator:");
            }

            virtual void value()
            {
                CONTEXT("When evaluating TransferOperator:");
                //TODO: think of new way to encode BCs -> separate operator?
                DenseVector<unsigned long> dummy(1ul);
                //std::cout << _left << _right << std::endl;
                TransferType_::value(_left, _right, dummy, _mat);
                //std::cout << _left << _right << std::endl;
            }

            virtual ~TransferOperator()
            {
                CONTEXT("When destroying TransferOperator:");
            }

        private:
            VectorType_ _left;
            VectorType_ _right;
            MatrixType_ _mat;
    };
}
#endif
