/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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


#include <honei/math/superlu.hh>
#include <honei/util/unittest.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/util/stringify.hh>
#include <honei/math/conjugate_gradients.hh>
#include <iostream>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class SuperLUTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _ref_f;
    public:
        SuperLUTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string ref_file) :
            BaseTest("Super LU Test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _ref_f = ref_file;
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _m_f;

            SparseMatrixELL<DT1_> smatrix(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT1_(0)));

            std::string filename_3(HONEI_SOURCEDIR);
            filename_3 += "/honei/math/testdata/";
            filename_3 += _ref_f;
            DenseVector<DT1_> ref(VectorIO<io_formats::EXP>::read_vector(filename_3, DT1_(0)));


            DenseVector<DT1_> result(rhs.size(), 4711);
            SuperLU::value(smatrix, rhs, result);

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
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref[i], eps);
                //std::cout<<result[i]<<"  "<<result_result[i]<<std::endl;
            }
        }
};

SuperLUTestSparseELL<tags::CPU, double> superlu_test_sparse_ell_double_1("double", "l2/area51_full_0.ell", "l2/area51_rhs_0", "l2/area51_sol_0");
SuperLUTestSparseELL<tags::CPU, double> superlu_test_sparse_ell_double_2("double", "poisson_advanced/q2_sort_0/A_7.ell", "poisson_advanced/q2_sort_0/rhs_7", "poisson_advanced/q2_sort_0/sol_7");
