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
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }
            SparseMatrixELL<DT_> smatrix(matrix);
            DenseVector<DT_> y1(Product<Tag_>::value(smatrix, x));
            DenseVector<DT_> y2(Product<tags::CPU>::value(matrix, x));

            y1.lock(lm_read_only);
            TEST_CHECK_EQUAL(y1, y2);
            y1.unlock(lm_read_only);
        }
};
MMIOTest<float, tags::CPU> mmio_test_float_dense("float");
MMIOTest<float, tags::GPU::CUDA> cuda_mmio_test_float_dense("float");

