/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/sparse_matrix.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>
#include <honei/mpi/sparse_matrix_csr_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/util/unittest.hh>
#include <honei/math/matrix_io.hh>
#include <honei/mpi/matrix_io_mpi.hh>

#include <string>
#include <limits>
#include <cmath>
#include <iostream>


using namespace honei;
using namespace tests;


template <typename DataType_>
class SparseMatrixMPIIOQuickTest :
    public QuickTest
{
    public:
        SparseMatrixMPIIOQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_mpi_io_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/poisson_advanced4/sort_0/A_4.ell";
            SparseMatrixELL<DataType_> saell(MatrixIO<io_formats::ELL>::read_matrix(filename, DataType_(0)));
            SparseMatrixELLMPI<DataType_> ref_aell(saell);
            SparseMatrixELLMPI<DataType_> aell(MatrixIOMPI_ELL<io_formats::ELL>::read_matrix(filename, DataType_(0)));

            TEST_CHECK_EQUAL(aell, ref_aell);

            SparseMatrixCSR<DataType_> sacsr(MatrixIO<io_formats::ELL>::read_matrix(filename, DataType_(0)));
            SparseMatrixCSRMPI<DataType_> ref_acsr(sacsr);
            SparseMatrixCSRMPI<DataType_> acsr(MatrixIOMPI_CSR<io_formats::ELL>::read_matrix(filename, DataType_(0)));

            TEST_CHECK_EQUAL(acsr, ref_acsr);
        }
};
SparseMatrixMPIIOQuickTest<double> sparse_matrix_mpi_io_quick_test_double("double");
