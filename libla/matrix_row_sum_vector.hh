/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH
#define LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH 1

#include <libla/dense_vector.hh>
#include <libla/matrix.hh>
#include <libla/vector_element_sum.hh>

#include <tr1/memory>

namespace pg512
{
    /**
     * A MatrixRowSumVector yields a Vector of which each element holds the sum
     * of the matrix's corresponding row's elements.
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct MatrixRowSumVector
    {
        static std::tr1::shared_ptr<DenseVector<DataType_> > value(const RowAccessMatrix<DataType_> & matrix)
        {
            std::tr1::shared_ptr<DenseVector<DataType_> > result(new DenseVector<DataType_>(matrix.rows()));

            for (typename Vector<DataType_>::ElementIterator i(result->begin_elements()), i_end(result->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = VectorElementSum<DataType_, Tag_>::value(matrix[i.index()]);
            }

            return result;
        }
    };
}

#endif
