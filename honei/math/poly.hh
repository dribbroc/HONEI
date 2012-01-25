/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#pragma once
#ifndef MATH_GUARD_POLY_HH
#define MATH_GUARD_POLY_HH 1

#include<honei/la/sparse_matrix.hh>
#include<honei/la/difference.hh>
#include<honei/la/sum.hh>
#include<honei/la/product.hh>
#include<honei/la/scale.hh>

namespace honei
{
    struct Poly
    {
        template<typename DT_>
        static inline SparseMatrix<DT_> value(SparseMatrix<DT_>& A, unsigned long m, DT_ damp = DT_(1))
        {
            //construct identity matrix
            SparseMatrix<DT_> I(A.rows(), A.columns());
            //construct D for diagonally scaled A
            SparseMatrix<DT_> D(A.rows(), A.columns());

            SparseMatrix<DT_> result(A.rows(), A.columns());

            for(unsigned long i(0) ; i < A.rows() ; ++i)
            {
                I(i, i, DT_(1));
                D(i, i, A(i, i) > std::numeric_limits<DT_>::epsilon() ? DT_(1) / A(i, i) : DT_(1) / std::numeric_limits<DT_>::epsilon());
            }
            SparseMatrix<DT_> DA(Product<tags::CPU>::value(D, A));

            SparseMatrix<DT_> I_minus_DA(A.rows(), A.columns());
            Difference<tags::CPU>::value(I_minus_DA, I, DA);

            //perform:

            SparseMatrix<DT_> temp(A.rows(), A.columns());
            //i := 0 => (I-DA)^0 = I
            temp = I.copy();
            result = temp.copy();

            //i > 0
            for(unsigned long k(1) ; k < m ; ++k)
            {
                if (k == 1)
                    temp = I_minus_DA.copy();
                else
                    temp = Product<tags::CPU>::value(temp.copy(), I_minus_DA);

                Sum<tags::CPU>::value(result, temp);
            }
            //multiply with D
            result = Product<tags::CPU>::value(result.copy(), D);
            return result;
        }
    };
}

#endif
