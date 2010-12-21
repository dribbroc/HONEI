/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HONEI_DUNE
#include <honei/math/conjugate_gradients.hh>
#include <honei/math/jacobi.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#ifdef HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif
#include <dune/istl/io.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace Dune;

template <typename Tag_, typename DT1_>
class DuneRegressionTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _r_f, _i_f;
    public:
        DuneRegressionTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string res_file,
                std::string init_file) :
            BaseTest("Dune regression test (sparse ELL system)<" + tag + ">")
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
            filename += "/honei/math/testdata/";
            filename += _m_f;
            SparseMatrix<DT1_> tsmatrix(MatrixIO<io_formats::M, SparseMatrix<DT1_> >::read_matrix(filename));
            SparseMatrixELL<DT1_> smatrix(tsmatrix);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));

            DenseVector<DT1_> diag_inverted(tsmatrix.rows(), DT1_(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = DT1_(1)/tsmatrix(i, i);
            }
            SparseMatrix<DT1_> tdifference(tsmatrix.copy());
            for(unsigned long i(0) ; i < tsmatrix.rows() ; ++i)
            {
                     tdifference(i, i) = DT1_(0);
            }
            SparseMatrixELL<DT1_> difference(tdifference);

            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _i_f;
            DenseVector<DT1_> init(VectorIO<io_formats::EXP>::read_vector(filename_4, DT1_(0)));
            DenseVector<DT1_> result(init.copy());

            unsigned long used_iters(0);
            ConjugateGradients<Tag_, JAC>::value(smatrix, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));
            //Jacobi<Tag_>::value(smatrix, difference, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));


            // DUNE

            typedef FieldMatrix<DT1_, 1, 1> M;
            typedef FieldVector<DT1_, 1> V;
            typedef BCRSMatrix<M> SM;
            typedef BlockVector<V> BV;

            SparseMatrix<DT1_> smatrix2(smatrix);
            SM smatrix_dune(smatrix2.rows(), smatrix2.columns(), SM::random);
            for (unsigned long i(0) ; i < smatrix2.rows() ; i++)
            {
                smatrix_dune.setrowsize(i, smatrix2[i].used_elements());
            }
            smatrix_dune.endrowsizes();
            for (typename SparseMatrix<DT1_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                smatrix_dune.addindex(i.row(), i.column());
            }
            smatrix_dune.endindices();
            for (typename SparseMatrix<DT1_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
            {
                smatrix_dune[i.row()][i.column()] = *i;
            }
            BV result_dune(result.size(), 0);
            BV rhs_dune(rhs.size(),0);
            for (unsigned long i(0) ; i < rhs.size() ; ++i)
            {
                rhs_dune[i] = rhs[i];
                result_dune[i] = init[i];
            }
            InverseOperatorResult irs;

#ifdef HAVE_SUPERLU
            SuperLU<SM> slu(smatrix_dune, true);
            slu.apply(result_dune, rhs_dune, irs);
#else
            typedef SeqJac<SM, BV, BV> PREC;
            PREC prec(smatrix_dune, 1, 1);
            MatrixAdapter<SM,BV,BV> op(smatrix_dune);
            CGSolver<BV> cg(op,prec,1E-13, 1000, 2);
            cg.apply(result_dune, rhs_dune, irs);
#endif

            result.lock(lm_read_only);
            result.unlock(lm_read_only);
            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);
            for (unsigned long i(0) ; i < result.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result_dune[i], result[i], eps);
                //std::cout<<result[i]<<"  "<<result_dune[i]<<std::endl;
            }

            //FEATFLOW
            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/";
            filename_3 += _r_f;
            DenseVector<DT1_> result_feat(VectorIO<io_formats::EXP>::read_vector(filename_3, DT1_(0)));


            std::cout << "Comparing with FEATFLOW2: eps= " << eps << std::endl;
            for(unsigned long i(0) ; i < result.size() ; ++i)
            {
                if(fabs(result[i] - result_feat[i]) > eps)
                    std::cout << std::setprecision(11) << result[i] << " " << result_feat[i] << " at index " << i << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], result_feat[i], eps);
            }
        }
};

#ifdef HONEI_SSE
DuneRegressionTestSparseELL<tags::CPU::SSE, double> sse_dune_regression_test_double_sparse_ell2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_sol_0", "l2/area51_init_0");
//DuneRegressionTestSparseELL<tags::CPU::SSE, double> sse_dune_regression_test_double_sparse_ell8("double", "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_sol_0", "l8/area51_init_0");
DuneRegressionTestSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_dune_regression_test_double_sparse_ell2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_sol_0", "l2/area51_init_0");
//DuneRegressionTestSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_dune_regression_test_double_sparse_ell8("double", "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_sol_0", "l8/area51_init_0");
#endif
#ifdef HONEI_CUDA_DOUBLE
DuneRegressionTestSparseELL<tags::GPU::CUDA, double> cuda_dune_regression_test_double_sparse_ell2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_sol_0", "l2/area51_init_0");
//DuneRegressionTestSparseELL<tags::GPU::CUDA, double> cuda_dune_regression_test_double_sparse_ell8("double", "l8/area51_full_0.m", "l8/area51_rhs_0", "l8/area51_sol_0", "l8/area51_init_0");
#endif

#endif
