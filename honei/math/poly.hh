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
        static inline SparseMatrix<DT_> value(SparseMatrix<DT_>& A, unsigned long m, DT_ damp)
        {
            //construct identity matrix
            SparseMatrix<DT_> A_damped(A.copy());
            Scale<tags::CPU>::value(A_damped, damp);
            SparseMatrix<DT_> I(A.rows(), A.columns());
            for(unsigned long i(0) ; i < A.rows() ; ++i)
            {
                I(i,i,DT_(1));
            }

            SparseMatrix<DT_> I_minus_A(A.rows(), A.columns());
            SparseMatrix<DT_> result(A.rows(), A.columns());
            Difference<tags::CPU>::value(I_minus_A, I, A_damped);

            //perform
            for(unsigned long k(0) ; k < m ; ++k)
            {

                if(k == 0)
                    Sum<tags::CPU>::value(result, I);
                else
                {
                    SparseMatrix<DT_> temp(I_minus_A.copy());
                    for(unsigned long i(0) ; i < k ; ++i)
                    {
                        temp = Product<tags::CPU>::value(temp, I_minus_A);
                    }

                    Sum<tags::CPU>::value(result, temp);
                }
            }
            return result;
        }
    };
}

#endif
