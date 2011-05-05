/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007-2011 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/bi_conjugate_gradients_stabilised.hh>
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
class PBiCGStabELLTEST:
    public BaseTest
{
    private:
        std::string _m_f, _precon_f, _v_f, _r_f, _i_f;
    public:
        PBiCGStabELLTEST(const std::string & tag,
                std::string m_file,
                std::string p_file,
                std::string v_file,
                std::string res_file,
                std::string init_file) :
            BaseTest("Preconditioned BiCGStab solver test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _precon_f = p_file;
            _v_f = v_file;
            _r_f = res_file;
            _i_f = init_file;
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _m_f;
            SparseMatrixELL<DT1_> smatrix2(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));

            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/";
            filename_3 += _precon_f;
            SparseMatrixELL<DT1_> smatrix3(MatrixIO<io_formats::ELL>::read_matrix(filename_3, DT1_(0)));

            /*SparseMatrix<DT1_> tsmatrix3(smatrix2.rows(), smatrix2.columns());
            for(unsigned long i(0) ; i < smatrix2.rows() ; ++i)
            {
                tsmatrix3(i , i) = 0.7 * 1/smatrix2(i, i);
            }
            SparseMatrixELL<DT1_> smatrix3(tsmatrix3);
*/
            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _i_f;
            DenseVector<DT1_> result(VectorIO<io_formats::EXP>::read_vector(filename_4, DT1_(0)));

            unsigned long used_iters(0);
            PBiCGStab<Tag_, methods::VAR>::value(smatrix2, rhs, result, smatrix3, 1000ul, used_iters, DT1_(1e-14));

            std::string filename_5(HONEI_SOURCEDIR);
            filename_5 += "/honei/math/testdata/";
            filename_5 += _r_f;
            DenseVector<DT1_> ref_result(VectorIO<io_formats::EXP>::read_vector(filename_5, DT1_(0)));

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
PBiCGStabELLTEST<tags::CPU, double> pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
//PBiCGStabELLTEST<tags::CPU::MultiCore, double> mc_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
#ifdef HONEI_SSE
PBiCGStabELLTEST<tags::CPU::SSE, double> sse_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
PBiCGStabELLTEST<tags::CPU::MultiCore::SSE, double> mcsse_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
PBiCGStabELLTEST<tags::GPU::CUDA, double> cuda_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
#endif
#endif
#ifdef HONEI_OPENCL
PBiCGStabELLTEST<tags::OpenCL::CPU, double> ocl_cpu_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
#ifdef HONEI_CUDA_DOUBLE
PBiCGStabELLTEST<tags::OpenCL::GPU, double> ocl_gpu_pbicgstab_test_double_sparse_ell("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/A_7_spai.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7", "poisson_advanced/sort_0/init_7");
#endif
#endif
