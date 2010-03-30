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
#include <unittest/unittest.hh>
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
        std::string _m_f, _v_f, _r_f;
    public:
        DuneRegressionTestSparseELL(const std::string & tag, std::string m_file, std::string v_file, std::string res_file) :
            BaseTest("Dune regression test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _r_f = res_file;
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _m_f;
            unsigned long non_zeros(MatrixIO<io_formats::M>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DT1_> data(non_zeros);

            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrix<DT1_> tsmatrix(rows, columns, r, c, data);
            SparseMatrixELL<DT1_> smatrix(tsmatrix);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(rows, DT1_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            DenseVector<DT1_> diag_inverted(rows, DT1_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DT1_(1)/data[i];
            }
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                     data[i] = DT1_(0);
            }
            SparseMatrix<DT1_> tdifference(rows, columns, r, c, data);
            SparseMatrixELL<DT1_> difference(tdifference);

            DenseVector<DT1_> result(rhs.size(), DT1_(0));
            DenseVector<DT1_> result_c(result.copy());
            Defect<Tag_>::value(result, rhs, smatrix, result_c);
            ConjugateGradients<Tag_, JAC>::value(smatrix, rhs, result, diag_inverted, 10000ul);
            //Jacobi<Tag_>::value(smatrix, difference, rhs, result, diag_inverted, 10000ul);


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
            result_dune=0;
            BV rhs_dune(rhs.size(),0);
            for (unsigned long i(0) ; i < rhs.size() ; ++i)
            {
                rhs_dune[i] = rhs[i];
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
            for (unsigned long i(0) ; i < result.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result_dune[i], result[i], 1e-8);
                //std::cout<<result[i]<<"  "<<result_dune[i]<<std::endl;
            }

            //FEATFLOW
            DenseVector<DT1_> result_feat(rows, DT1_(0));
            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/";
            filename_3 += _r_f;
            VectorIO<io_formats::EXP>::read_vector(filename_3, result_feat);
            std::cout << "Comparing with FEATFLOW2: " << std::endl;
            for(unsigned long i(0) ; i < result.size() ; ++i)
            {
                //std::cout<<result[i]<<"..."<<result_dune[i]<<"..."<<result_feat[i]<<std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result_dune[i], result_feat[i], 1e-4);
            }
        }
};

DuneRegressionTestSparseELL<tags::CPU::SSE, double> sse_dune_regression_test_double_sparse_ell2("double", "l2/poisson_full.m", "l2/poisson_rhs", "l2/poisson_sol");
//DuneRegressionTestSparseELL<tags::CPU::SSE, double> sse_dune_regression_test_double_sparse_ell8("double", "l8/poisson_full.m", "l8/poisson_rhs", "l8/poisson_sol");
DuneRegressionTestSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_dune_regression_test_double_sparse_ell2("double", "l2/poisson_full.m", "l2/poisson_rhs", "l2/poisson_sol");
DuneRegressionTestSparseELL<tags::CPU::MultiCore::SSE, double> mc_sse_dune_regression_test_double_sparse_ell8("double", "l8/poisson_full.m", "l8/poisson_rhs", "l8/poisson_sol");
#ifdef HONEI_CUDA_DOUBLE
DuneRegressionTestSparseELL<tags::GPU::CUDA, double> cuda_dune_regression_test_double_sparse_ell2("double", "l2/poisson_full.m", "l2/poisson_rhs", "l2/poisson_sol");
DuneRegressionTestSparseELL<tags::GPU::CUDA, double> cuda_dune_regression_test_double_sparse_ell8("double", "l8/poisson_full.m", "l8/poisson_rhs", "l8/poisson_sol");
#endif

#endif
