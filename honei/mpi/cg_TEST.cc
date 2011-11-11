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

#include <honei/la/dense_vector.hh>
#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/util/unittest.hh>
#include <honei/math/cg.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/time_stamp.hh>

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>


using namespace honei;
using namespace tests;


template <typename Tag_, typename DT1_>
class CGSolverTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _r_f, _i_f;
    public:
        CGSolverTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string res_file,
                std::string init_file) :
            BaseTest("MPI CG solver test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _r_f = res_file;
            _i_f = init_file;
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/poisson_advanced2/q2_sort_0/";
            filename += _m_f;
            SparseMatrixELL<DT1_> smatrix2(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));
            SparseMatrix<DT1_> ssmatrix2(smatrix2);
            SparseMatrixELLMPI<DT1_> matrix2(ssmatrix2);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/poisson_advanced2/q2_sort_0/";
            filename_2 += _v_f;
            DenseVector<DT1_> srhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));
            DenseVectorMPI<DT1_> rhs(srhs);

            DenseVector<DT1_> sdiag_inverted(smatrix2.rows(), DT1_(0));
            for(unsigned long i(0) ; i < sdiag_inverted.size() ; ++i)
            {
                    sdiag_inverted[i] = DT1_(0.7)/smatrix2(i, i);
            }
            DenseVectorMPI<DT1_> diag_inverted(sdiag_inverted);

            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/poisson_advanced2/q2_sort_0/";
            filename_4 += _i_f;
            DenseVector<DT1_> sresult(VectorIO<io_formats::EXP>::read_vector(filename_4, DT1_(0)));
            DenseVectorMPI<DT1_> result(sresult);

            unsigned long used_iters(4711);
            TimeStamp at, bt;
            at.take();
            CG<Tag_, methods::VAR>::value(matrix2, diag_inverted, rhs, result, 1000ul, used_iters, 1e-6);
            bt.take();
            if (mpi::mpi_comm_rank() == 0)
            {
                std::cout<<"Used iters: "<<used_iters<<std::endl;
                std::cout<<"TOE: "<<bt.total()-at.total()<<std::endl;
            }

            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/poisson_advanced2/q2_sort_0/";
            filename_3 += _r_f;
            DenseVector<DT1_> sref_result(VectorIO<io_formats::EXP>::read_vector(filename_3, DT1_(0)));
            DenseVectorMPI<DT1_> ref_result(sref_result);



            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);

            if (mpi::mpi_comm_rank() == 0)
            {
                std::cout << "Comparing with FEATFLOW2: eps= " << eps << std::endl;
            }
            for(unsigned long i(0) ; i < result.local_size() ; ++i)
            {
                if(fabs(result[i] - ref_result[i]) > eps)
                    std::cout << std::setprecision(11) << result[i] << " " << ref_result[i] << " at index " << i << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], eps);
            }
        }
};
#ifdef HONEI_SSE
CGSolverTestSparseELL<tags::CPU::SSE, double> cgs_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#else
CGSolverTestSparseELL<tags::CPU, double> cgs_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#endif
