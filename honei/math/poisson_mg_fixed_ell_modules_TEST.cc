/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
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

#define SOLVER_VERBOSE_L2
#define SOLVER_VERBOSE

#include <honei/math/multigrid.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class PoissonTestMGModulesSparseELL:
    public BaseTest
{
    private:
        unsigned long _size;
        std::string _res_f;

    public:
        PoissonTestMGModulesSparseELL(const std::string & tag) :
            BaseTest("Poisson test for itrerative LES solvers , MG Modules (ELLPACK vs. Q1 system)<" + tag + ">")
    {
        register_tag(Tag_::name);
    }
        virtual void run() const
        {
            for (unsigned long i(0) ; i <= 6; ++i)
            {
                unsigned long size((unsigned long)(((unsigned long)pow((DT1_)2, (DT1_)i) + 1) * ((unsigned long)pow((DT1_)2, (DT1_)i) + 1)));
                if(i == 0)
                    size = 9;

                std::cout << "Level: " << size << std::endl;

                DenseVector<DT1_> rhs(size, DT1_(0));
                DenseVector<DT1_> dummy_band(size, DT1_(0));


                FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(rhs);

                BandedMatrixQ1<DT1_> banded_matrix(size, dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy(), dummy_band.copy());
                FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(banded_matrix);
                SparseMatrix<DT1_> sm(banded_matrix);
                SparseMatrixELL<DT1_> sparse_matrix(sm);

                DenseVector<DT1_> banded_result(size, DT1_(0));
                DenseVector<DT1_> sparse_result(size, DT1_(0));

                DenseVector<DT1_> diag(banded_matrix.band(DD).copy());

                Jacobi<Tag_>::value(banded_matrix, rhs, banded_result, DT1_(0.7), diag);
                Jacobi<Tag_>::value(sparse_matrix, rhs, sparse_result, DT1_(0.7), diag);

                banded_result.lock(lm_read_only);
                banded_result.unlock(lm_read_only);

                sparse_result.lock(lm_read_only);
                sparse_result.unlock(lm_read_only);
                TEST_CHECK_EQUAL(sparse_result, banded_result);

                //----------------------------------------------------------------------
                DenseVector<DT1_> banded_result_3(size, DT1_(0));
                DenseVector<DT1_> sparse_result_3(size, DT1_(0));
                ConjugateGradients<Tag_, NONE>::value(banded_matrix, rhs, banded_result_3, std::numeric_limits<DT1_>::epsilon());
                ConjugateGradients<Tag_, NONE>::value(sparse_matrix, rhs, sparse_result_3, std::numeric_limits<DT1_>::epsilon());
                banded_result_3.lock(lm_read_only);
                banded_result_3.unlock(lm_read_only);

                sparse_result_3.lock(lm_read_only);
                sparse_result_3.unlock(lm_read_only);
                for(unsigned long j(0) ; j < banded_result_3.size() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS(sparse_result_3[j], banded_result_3[j], std::numeric_limits<DT1_>::epsilon() * DT1_(1e3));
                //----------------------------------------------------------------------

                DenseVector<DT1_> banded_result_2(size, DT1_(0));
                DenseVector<DT1_> sparse_result_2(size, DT1_(0));
                Jacobi<Tag_>::value(banded_result_2, banded_matrix, rhs, banded_result_2, 1, DT1_(0.7), diag);
                Jacobi<Tag_>::value(sparse_result_2, sparse_matrix, rhs, sparse_result_2, 1, DT1_(0.7), diag);
                banded_result_2.lock(lm_read_only);
                banded_result_2.unlock(lm_read_only);

                sparse_result_2.lock(lm_read_only);
                sparse_result_2.unlock(lm_read_only);
                for(unsigned long j(0) ; j < banded_result_2.size() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS(sparse_result_2[j], banded_result_2[j], std::numeric_limits<DT1_>::epsilon());
                //----------------------------------------------------------------------

            }
        }
};
PoissonTestMGModulesSparseELL<tags::CPU, float> poisson_test_mg_banded_float("float");
PoissonTestMGModulesSparseELL<tags::CPU, double> poisson_test_mg_banded_double("double");
#ifdef HONEI_SSE
PoissonTestMGModulesSparseELL<tags::CPU::SSE, float> poisson_test_mg_banded_float_sse("SSE float");
PoissonTestMGModulesSparseELL<tags::CPU::SSE, double> poisson_test_mg_banded_double_sse("SSE double");
#endif
#ifdef HONEI_CUDA
PoissonTestMGModulesSparseELL<tags::GPU::CUDA, float> poisson_test_mg_banded_float_cuda("CUDA float");
#ifdef HONEI_CUDA_DOUBLE
PoissonTestMGModulesSparseELL<tags::GPU::CUDA, double> poisson_test_mg_banded_double_cuda("CUDA double");
#endif
#endif
