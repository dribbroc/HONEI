/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <iostream>

#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/mpi/sparse_matrix_csr_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/util/unittest.hh>
#include <honei/math/matrix_io.hh>

#include <string>
#include <limits>
#include <cmath>


using namespace honei;
using namespace tests;


template <typename DT_>
class SparseMatrixCSRMPIQuickTest :
    public QuickTest
{
    public:
        SparseMatrixCSRMPIQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_csr_mpi_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::string dir(HONEI_SOURCEDIR);
            std::string file (dir + "/honei/math/testdata/poisson_advanced2/sort_0/");
            file += "prol_4";
            file += ".ell";
            SparseMatrixELL<DT_> aell(MatrixIO<io_formats::ELL>::read_matrix(file, DT_(0)));
            SparseMatrixCSR<DT_> acsr(aell);

            SparseMatrix<DT_> as(acsr);
            //SparseMatrixCSRMPI<DT_> a(as);
            SparseMatrixCSRMPI<DT_> a(acsr);


            SparseMatrixCSRMPI<DT_> sm0(as);
            SparseMatrixCSRMPI<DT_> sm1(sm0);
            SparseMatrixCSRMPI<DT_> sm2(sm1.copy());

            //TEST_CHECK_EQUAL(sm2, sm0);

            for (unsigned long i(0) ; i < sm2.local_rows() ; ++i)
                for (unsigned long j(0) ; j < sm2.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL(sm0(i, j), acsr(i + sm2.offset(), j));
                }
        }
};
SparseMatrixCSRMPIQuickTest<double> sparse_matrix_csr_mpi_quick_test_double("double");

template <typename DT_>
class SparseMatrixPartCSRMPIQuickTest :
    public QuickTest
{
    public:
        SparseMatrixPartCSRMPIQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_part_csr_mpi_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            SparseMatrix<DT_> as(12, 13);
            for (typename SparseMatrix<DT_>::ElementIterator i(as.begin_elements()) ; i < as.end_elements() ; ++i)
            {
                if (i.index() % 7 == 0)
                    as(i.row(), i.column(), DT_(i.index()) / DT_(mpi::mpi_comm_rank() + 1));
            }

            SparseMatrixCSRMPI<DT_> sm0(as, mpi::mpi_comm_size() * as.rows());
            SparseMatrixCSRMPI<DT_> sm1(sm0);
            SparseMatrixCSRMPI<DT_> sm2(sm1.copy());
            //TEST_CHECK_EQUAL(sm2, sm0);

            for (unsigned long i(0) ; i < sm2.local_rows() ; ++i)
                for (unsigned long j(0) ; j < sm2.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL(sm0(i, j), as(i, j));
                }
        }
};
SparseMatrixPartCSRMPIQuickTest<double> sparse_matrix_part_csr_mpi_quick_test_double("double");
