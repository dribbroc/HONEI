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

#ifndef MATH_GUARD_OPERATOR_HH
#define MATH_GUARD_OPERATOR_HH 1

#include <honei/math/defect.hh>
#include <honei/la/norm.hh>
#include <honei/la/sum.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/math/cg.hh>

namespace honei
{
    class Operator
    {
        public:
            virtual void value() = 0;
            virtual ~Operator()
            {
            }
    };

    class OperatorList
    {
        private:
            std::vector<Operator*> _ops;

        public:
            ~OperatorList()
            {
                for (unsigned long i(0) ; i < _ops.size() ; ++i)
                    delete _ops.at(i);
            }

            void push_back(Operator * op)
            {
                _ops.push_back(op);
            }

            void value()
            {
                for (unsigned long i(0) ; i < _ops.size() ; ++i)
                    _ops.at(i)->value();
            }

            Operator* operator[](unsigned long i)
            {
                return _ops.at(i);
            }
    };

    template<typename Tag_, typename MatType_, typename DT_>
    class DefectOperator : public Operator
    {
        public:
            DefectOperator(DenseVector<DT_>& y, DenseVector<DT_>& b, MatType_& A, DenseVector<DT_>& x) :
                _y(y),
                _b(b),
                _A(A),
                _x(x)
            {
            }

            virtual void value()
            {
                Defect<Tag_>::value(_y, _b, _A, _x);
            }

            virtual ~DefectOperator()
            {
            }

        private:
            DenseVector<DT_> _y;
            DenseVector<DT_> _b;
            MatType_ _A;
            DenseVector<DT_> _x;
    };

    template<typename Tag_, VectorNormType normtype, typename DT_, bool root>
    class NormOperator : public Operator
    {
        public:
            NormOperator(DT_& y, DenseVector<DT_>& x) :
                _x(x),
                _y(y)
        {
        }

            virtual void value()
            {
                _y = Norm<normtype, root, Tag_>::value(_x);
            }

            virtual ~NormOperator()
            {
            }

        private:
            DenseVector<DT_> _x;
            DT_& _y;
    };

    template<typename Tag_, typename DT_>
    class SumOperator : public Operator
    {
        public:
            SumOperator(DenseVector<DT_>& y, DenseVector<DT_>& x) :
                _x(x),
                _y(y)
        {
        }

            virtual void value()
            {
                Sum<Tag_>::value(_y, _x);
            }

            virtual ~SumOperator()
            {
            }

        private:
            DenseVector<DT_> _x;
            DenseVector<DT_> _y;
    };

    template<typename Tag_, typename DT_>
    class ScaledSumOperator : public Operator
    {
        public:
            ScaledSumOperator(DenseVector<DT_>& y, DenseVector<DT_>& x, DT_ alpha) :
                _x(x),
                _y(y),
                _alpha(alpha)
        {
        }

            virtual void value()
            {
                ScaledSum<Tag_>::value(_y, _x, _alpha);
            }

            virtual ~ScaledSumOperator()
            {
            }

        private:
            DenseVector<DT_> _x;
            DenseVector<DT_> _y;
            DT_ _alpha;
    };

    template<typename Tag_, typename SolverType_, typename MatrixType_, typename DT_>
    class SolverOperator : public Operator
    {
        //TODO: how to build in preconditioning -> standard-value?
        public:
            SolverOperator(MatrixType_ & A,
                           DenseVector<DT_> & b,
                           DenseVector<DT_> & x,
                           unsigned long max_iters,
                           unsigned long & used_iters,
                           DT_ eps_relative = 1e-8) :
                _A(A),
                _b(b),
                _x(x),
                _max_iters(max_iters),
                _used_iters(used_iters),
                _eps_relative(eps_relative)
            {
            }

            virtual void value()
            {
                SolverType_::value(_A, _b, _x, _max_iters, _used_iters, _eps_relative);
            }

        private:
            MatrixType_ _A;
            DenseVector<DT_> _b;
            DenseVector<DT_> _x;
            unsigned long _max_iters;
            unsigned long & _used_iters;
            DT_ _eps_relative;
    };
}
#endif
