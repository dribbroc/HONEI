/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH
#define LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH 1

#include <libla/matrix.hh>
#include <libla/vector_sum.hh>

#include <tr1/memory>

namespace pg512
{
    /**
     * A MatrixRowSumVector yields a column-vector of which each element holds the sum
     * of the matrix's corresponding row-vector.
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct MatrixRowSumVector
    {
        static std::tr1::shared_ptr<DenseVector<DataType_> > value(const Matrix<DataType_> & matrix)
        {
            std::tr1::shared_ptr<DenseVector<DataType_> > result(new DenseVector<DataType_>(matrix.rows()));

            for (unsigned long r(0) ; r < matrix.rows() ; ++r)
            {
                (*result)[r] = VectorSum<DataType_, Tag_>::value(matrix[r]);
            }

            return result;
        }
    };
}

#endif
