/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#define SOLVER_VERBOSE_L2

#include <honei/math/conjugate_gradients.hh>
#include <honei/math/jacobi.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/math/dune_solvers.hh>
#include <iostream>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class DuneConjugateGradientsTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _i_f;
    public:
        DuneConjugateGradientsTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string init_file) :
            BaseTest("Dune conjugate gradients test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _i_f = init_file;
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

            DenseVector<DT1_> init(rhs.size(), DT1_(0));
            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _i_f;
            VectorIO<io_formats::EXP>::read_vector(filename_4, init);
            DenseVector<DT1_> result(init.copy());

            unsigned long used_iters(0);
            ConjugateGradients<Tag_, JAC>::value(smatrix, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));
            //Jacobi<Tag_>::value(smatrix, difference, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));

            // DUNE

            DenseVector<DT1_> dune_result(init.copy());
            DuneConjugateGradients<JAC>::value(smatrix, rhs, dune_result, diag_inverted, 10000ul, DT1_(1e-8));

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
                TEST_CHECK_EQUAL_WITHIN_EPS(dune_result[i], result[i], eps);
                //std::cout<<result[i]<<"  "<<result_result[i]<<std::endl;
            }
        }
};

#ifdef HONEI_SSE
DuneConjugateGradientsTestSparseELL<tags::CPU::SSE, double> sse_dune_cg_test_sparse_ell_double_2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
DuneConjugateGradientsTestSparseELL<tags::GPU::CUDA, double> cuda_dune_cg_test_sparse_ell_double_2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
#endif
#endif

template <typename Tag_, typename DT1_>
class DuneLUTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _i_f;
    public:
        DuneLUTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string init_file) :
            BaseTest("Dune LU Test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _i_f = init_file;
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

            DenseVector<DT1_> init(rhs.size(), DT1_(0));
            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _i_f;
            VectorIO<io_formats::EXP>::read_vector(filename_4, init);
            DenseVector<DT1_> result(init.copy());

            unsigned long used_iters(0);
            ConjugateGradients<Tag_, JAC>::value(smatrix, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));
            //Jacobi<Tag_>::value(smatrix, difference, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));


            // DUNE

            DenseVector<DT1_> dune_result(init.copy());
            DuneLU::value(smatrix, rhs, dune_result);

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
                TEST_CHECK_EQUAL_WITHIN_EPS(dune_result[i], result[i], eps);
                //std::cout<<result[i]<<"  "<<result_result[i]<<std::endl;
            }
        }
};

#ifdef HONEI_SSE
DuneLUTestSparseELL<tags::CPU::SSE, double> sse_dune_lu_test_sparse_ell_double_2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
DuneLUTestSparseELL<tags::GPU::CUDA, double> cuda_dune_lu_test_sparse_ell_double_2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");
#endif
#endif

#endif
