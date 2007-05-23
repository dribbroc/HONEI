/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <libla/dense_matrix.hh>
#include <libutil/exception.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace pg512;

class DenseMatrixCreationTest :
    public BaseTest
{
    public:
        DenseMatrixCreationTest() :
            BaseTest("dense_matrix_creation_test")
        {
        }

        virtual void run() const
        {
            for (unsigned long size ; size < (1 << 10) ; size << 1)
            {
                std::tr1::shared_ptr<DenseMatrix<float> > dm(new DenseMatrix<float>(size, size, 0));
            }
        }
} dense_matrix_creation_test;
