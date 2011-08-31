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

#include <honei/math/cg.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <iomanip>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <honei/math/sainv.hh>
#include <honei/math/superlu.hh>
#include <honei/math/ri.hh>
#include <honei/math/bi_conjugate_gradients_stabilised.hh>
#include <honei/math/spai2.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class FeastTestSparseELLPrecon:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _r_f, _i_f;
    public:
        FeastTestSparseELLPrecon(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string res_file,
                std::string init_file) :
            BaseTest("Feast solver test (sparse ELL system)<" + tag + ">")
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
            filename += "/honei/math/testdata/feast_aniso/";
            filename += _m_f;
            SparseMatrixELL<DT1_> smatrix2(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/feast_aniso/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(VectorIO<io_formats::DV>::read_vector(filename_2, DT1_(0)));

            /*DenseVector<DT1_> diag_inverted(smatrix2.rows(), DT1_(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = DT1_(1.0)/smatrix2(i, i);
            }*/
            SparseMatrix<DT1_> sm(smatrix2);
            SparseMatrix<DT1_> Ms = SAINV<Tag_>::value(sm);
            /*SparseMatrix<DT1_> Ms(sm.copy());
            SPAI2<Tag_>::value(Ms, sm);*/
            SparseMatrixELL<DT1_> diag_inverted(Ms);

            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/feast_aniso/";
            filename_4 += _i_f;
            DenseVector<DT1_> result(VectorIO<io_formats::DV>::read_vector(filename_4, DT1_(0)));
            for (unsigned long i(0) ; i < result.size() ; ++i)
                if (result[i] != DT1_(1))
                    result[i] = DT1_(0);
            std::cout<<rhs;
            std::cout<<result;
            unsigned long used_iters;
            //CG<Tag_, methods::VAR>::value(smatrix2, diag_inverted, rhs, result, 10000ul, used_iters, DT1_(1e-8));
            PBiCGStab<Tag_, methods::VAR>::value(smatrix2, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));
            //CG<Tag_, methods::NONE>::value(smatrix2, rhs, rhs, result, 10000ul, used_iters, DT1_(1e-8));
            /*DenseVector<DT1_> temp_0(rhs.size());
            DenseVector<DT1_> temp_1(rhs.size());
            RISmoother<Tag_>::value(smatrix2, diag_inverted, rhs, result, temp_0, temp_1, 750ul); //200 */
            std::cout<<result;

            DenseVector<DT1_> ref_result(rhs.size());
            SuperLU::value(smatrix2, rhs, ref_result);

            /*std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/feast_aniso/";
            filename_3 += _r_f;
            DenseVector<DT1_> ref_result(VectorIO<io_formats::DV>::read_vector(filename_3, DT1_(0)));*/
            std::cout<<ref_result;

            result.lock(lm_read_only);
            //std::cout << result << std::endl;
            result.unlock(lm_read_only);

            std::cout << "Used iters: " << used_iters << std::endl;

            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);
            //eps *= DT1_(100);

            std::cout << "Comparing with FEATFLOW2: eps= " << eps << std::endl;
            for(unsigned long i(0) ; i < result.size() ; ++i)
            {
                if(fabs(result[i] - ref_result[i]) > eps)
                    std::cout << std::setprecision(11) << result[i] << " " << ref_result[i] << " at index " << i << std::endl;
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], eps);
            }
        }
};
FeastTestSparseELLPrecon<tags::CPU, double> feast_precon_test_double_sparse_ell("double", "A_4.ell", "289_rhs.dv", "289_sol.dv", "289_rhs.dv");
