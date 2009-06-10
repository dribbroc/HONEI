/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
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

#include <honei/math/matrix_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/la/product.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template<typename DT_, typename Tag_>
class MMIOTest:
    public BaseTest
{
    public:
        MMIOTest(const std::string & tag) :
            BaseTest("Matrix Market Format I/O test")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/5pt_10x10.mtx";
            unsigned long non_zeros(0);
            DenseMatrix<DT_> matrix = MatrixIO::read_matrix(filename, DT_(0), non_zeros);
            std::cout << matrix;
            std::cout << "Non zero elements: " << non_zeros << std::endl;

            DenseVector<DT_> x(matrix.rows());
            DenseVector<DT_> y(matrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }
            SparseMatrix<DT_> ssmatrix(matrix);
            SparseMatrixELL<DT_> smatrix(ssmatrix);
            DenseVector<DT_> y1(Product<Tag_>::value(y, smatrix, x));
            DenseVector<DT_> y2(Product<tags::CPU>::value(matrix, x));

            y1.lock(lm_read_only);
            TEST_CHECK_EQUAL(y1, y2);
            y1.unlock(lm_read_only);

            ///Test second read-in method:
            std::cout << "Testing sparse data read-in:" << std::endl;
            unsigned long non_zeros_2(MatrixIO::get_non_zeros(filename));
            TEST_CHECK_EQUAL(non_zeros, non_zeros_2);
            DenseVector<unsigned long> r(non_zeros_2);
            DenseVector<unsigned long> c(non_zeros_2);
            DenseVector<DT_> data(non_zeros_2);

            MatrixIO::read_matrix(filename, r, c, data);
            std::cout << "Row Indices:" << std::endl;
            std::cout << r;
            std::cout << "Column Indices:" << std::endl;
            std::cout << c;
            std::cout << "Data:" << std::endl;
            std::cout << data;

            for(unsigned long i(0) ; i < non_zeros_2 ; ++i)
            {
                TEST_CHECK_EQUAL(matrix[r[i]][c[i]] , data[i]);
            }

            unsigned long rows, columns, ax, bx;
            MatrixIO::get_sizes(filename, columns, rows, ax, bx);
            SparseMatrixELL<DT_> smatrix2(rows, columns, r, c, data);

            TEST_CHECK_EQUAL(smatrix2, smatrix);
        }
};
MMIOTest<float, tags::CPU> mmio_test_float_dense("float");
MMIOTest<double, tags::CPU> mmio_test_double_dense("double");
#ifdef HONEI_CUDA
MMIOTest<float, tags::GPU::CUDA> cuda_mmio_test_float_dense("float");
#ifdef HONEI_CUDA_DOUBLE
MMIOTest<double, tags::GPU::CUDA> cuda_mmio_test_double_dense("double");
#endif
#endif

