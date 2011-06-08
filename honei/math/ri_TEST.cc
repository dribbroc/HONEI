/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. LibMath is free software;
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

#include <honei/math/ri.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <iomanip>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class RITestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _r_f, _i_f;
    public:
        RITestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string res_file,
                std::string init_file) :
            BaseTest("RI smoother test (sparse ELL system)<" + tag + ">")
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
            filename += "/honei/math/testdata/poisson_advanced/sort_0/";
            filename += _m_f;
            SparseMatrixELL<DT1_> smatrix2(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/poisson_advanced/sort_0/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));

            DenseVector<DT1_> diag_inverted(smatrix2.rows(), DT1_(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = DT1_(0.7)/smatrix2(i, i);
            }

            /*SparseMatrix<DT1_> bla(smatrix2.rows(), smatrix2.columns());
            for(unsigned long i(0) ; i < bla.rows() ; ++i)
                for(unsigned long j(0) ; j < bla.columns() ; ++ j)
                {
                    if (smatrix2(i,j) != DT1_(0))
                        bla(i,j) = smatrix2(i,j);
                }
            DenseMatrix<DT1_> dmatrix(bla);*/
            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/poisson_advanced/sort_0/";
            filename_4 += _i_f;
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(filename_4, DT1_(0)));

            DenseVector<DT1_> temp_0(rhs.size());
            DenseVector<DT1_> temp_1(rhs.size());
            RISmoother<Tag_>::value(smatrix2, diag_inverted, rhs, result, temp_0, temp_1, 1000ul);

            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/poisson_advanced/sort_0/";
            filename_3 += _r_f;
            DenseVector<DT1_> ref_result(VectorIO<io_formats::EXP>::read_vector(filename_3, DT1_(0)));

            result.lock(lm_read_only);
            //std::cout << result << std::endl;
            result.unlock(lm_read_only);


            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);

            std::cout << "Comparing with FEATFLOW2: eps= " << eps << std::endl;
            for(unsigned long i(0) ; i < result.size() ; ++i)
            {
                if(fabs(result[i] - ref_result[i]) > eps)
                    std::cout << std::setprecision(11) << result[i] << " " << ref_result[i] << " at index " << i << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], eps);
            }
        }
};
RITestSparseELL<tags::CPU, double> ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#ifdef HONEI_SSE
RITestSparseELL<tags::CPU::SSE, double> sse_ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
RITestSparseELL<tags::CPU::MultiCore::SSE, double> mcsse_ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
RITestSparseELL<tags::GPU::CUDA, double> cuda_ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#endif
#endif
#ifdef HONEI_OPENCL
RITestSparseELL<tags::OpenCL::CPU, double> ocl_cpu_ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#ifdef HONEI_CUDA_DOUBLE
RITestSparseELL<tags::OpenCL::GPU, double> ocl_gpu_ri_test_double_sparse_ell("double", "A_4.ell", "rhs_4", "sol_4", "init_4");
#endif
#endif
