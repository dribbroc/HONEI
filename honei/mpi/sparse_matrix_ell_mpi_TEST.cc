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
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/util/unittest.hh>

#include <string>
#include <limits>
#include <cmath>
#include <iostream>


using namespace honei;
using namespace tests;


template <typename DataType_>
class SparseMatrixELLMPIQuickTest :
    public QuickTest
{
    public:
        SparseMatrixELLMPIQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_ell_mpi_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            mpi::mpi_init();
            int rank;
            mpi::mpi_comm_rank(&rank);
            int comm_size;
            mpi::mpi_comm_size(&comm_size);

            SparseMatrix<DataType_> src(711, 713);
            for (unsigned long i(0) ; i < src.rows() ; ++i)
                for (unsigned long j(0) ; j < src.columns() ; ++j)
                    if (i % 7 == 0)
                        src(i, j, i+10);


            SparseMatrixELLMPI<DataType_> sm0(src, rank, comm_size);
            SparseMatrixELLMPI<DataType_> sm1(sm0);
            SparseMatrixELLMPI<DataType_> sm2(sm1.copy());
            std::cout<<sm2.rows()<<" "<<sm2.offset()<<std::endl;

            for (unsigned long i(0) ; i < sm2.rows() ; ++i)
                for (unsigned long j(0) ; j < sm2.columns() ; ++j)
                    TEST_CHECK_EQUAL(sm2(i, j), src(i + sm2.offset(), j));


            mpi::mpi_finalize();
        }
};
SparseMatrixELLMPIQuickTest<double> sparse_matrix_ell_mpi_quick_test_double("double");
