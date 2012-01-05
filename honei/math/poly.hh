/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
            SparseMatrix<DT_> A_damped(A.copy());
            Scale<tags::CPU>::value(A_damped, damp);
            //construct identity matrix
            SparseMatrix<DT_> I(A.rows(), A.columns());
            //construct D for diagonally scaled A
            SparseMatrix<DT_> D(A.rows(), A.columns());

            SparseMatrix<DT_> result(A.rows(), A.columns());

            for(unsigned long i(0) ; i < A.rows() ; ++i)
            {
                result(i, i, DT_(1));
                I(i, i, DT_(1));
                D(i, i, DT_(1) / A(i, i));
            }
            SparseMatrix<DT_> DA(Product<tags::CPU>::value(D, A_damped));

            SparseMatrix<DT_> I_minus_A(A.rows(), A.columns());
            Difference<tags::CPU>::value(I_minus_A, I, DA);

            //perform
            for(unsigned long k(1) ; k < m ; ++k)
            {

                SparseMatrix<DT_> temp(result.copy());
                if (k == 1)
                    temp = I_minus_A.copy();
                temp = Product<tags::CPU>::value(temp, I_minus_A);

                Sum<tags::CPU>::value(result, temp);
            }
            Sum<tags::CPU>::value(result, I);
            result = Product<tags::CPU>::value(result, D);
            return result;
        }
    };
}

#endif
