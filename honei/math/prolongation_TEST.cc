/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/prolongation.hh>
#include <honei/math/reordering.hh>
#include <honei/math/prolongation_matrix.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/la/dense_matrix.hh>
#include <cmath>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class ProlongationTest:
    public BaseTest
{
    public:
        ProlongationTest(const std::string & tag) :
            BaseTest("Prolongate test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (float level(0) ; level < 10 ; ++level)
            {
                unsigned long N_fine((unsigned long)pow((pow(2, level + 1) + 1), 2));
                //unsigned long width_fine = (unsigned long)sqrt((double)N_fine);
                unsigned long N_coarse((unsigned long)pow((pow(2, level) + 1), 2));
                //unsigned long width_coarse = (unsigned long)sqrt((double)N_coarse);

                DenseVector<DT1_> fine(N_fine, DT1_(4711));
                DenseVector<DT1_> fine_ref(N_fine, DT1_(42));
                DenseVector<DT1_> coarse(N_coarse, DT1_(1));
                /*for (unsigned long i(0) ; i < coarse.size() ; ++i)
                {
                    coarse[i] = DT1_(i % 1000);
                }*/
                DenseVector<unsigned long> mask(8);
                for(unsigned long i(0) ; i < 8 ; ++i)
                {
                    mask[i] = 2ul;
                }

                SparseMatrix<DT1_> sm(1, 1);
                SparseMatrixELL<DT1_> dummy(sm);
                Prolongation<Tag_, methods::NONE>::value(fine, coarse, mask, dummy);
                Prolongation<tags::CPU, methods::NONE>::value(fine_ref, coarse, mask, dummy);
                fine.lock(lm_read_only);
                fine_ref.lock(lm_read_only);
                fine.unlock(lm_read_only);
                fine_ref.unlock(lm_read_only);
                std::cout << "At level: " << level + 1 << std::endl;
                for(unsigned long i(0) ; i < fine.size() ; ++i)
                {
                    if (fine[i] != fine_ref[i])
                        std::cout << "Not matching: " << i << std::endl;
                }
                TEST_CHECK_EQUAL(fine, fine_ref);
            }
        }
};
#ifdef HONEI_CUDA
ProlongationTest<tags::GPU::CUDA, float> cuda_prolongate_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
ProlongationTest<tags::GPU::CUDA, double> cuda_prolongate_test_double("double");
#endif
#endif

template <typename Tag_, typename DT1_>
class ProlongationPROLMATTest:
    public BaseTest
{
    private:
        unsigned long _fine_size;
        unsigned long _coarse_size;
    public:
        ProlongationPROLMATTest(const std::string & tag, unsigned long fine, unsigned long coarse) :
            BaseTest("Prolongation PROLMAT test <" + tag + ">")
        {
            _fine_size = fine;
            _coarse_size = coarse;
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> defect_coarse(_coarse_size, DT1_(1.234));
            DenseVector<DT1_> defect_fine(_fine_size);

            DenseVector<unsigned long> fine_index(_fine_size);
            DenseVector<unsigned long> coarse_index(_coarse_size);
            DenseVector<unsigned long> mask(8);
            for(unsigned long i(0) ; i < 8 ; ++i)
            {
                mask[i] = 2ul;
            }

            Reordering<Tag_, methods::NATURAL>::value(fine_index, coarse_index);

            SparseMatrix<DT1_> sm(_fine_size, _coarse_size);
            ProlongationMatrix<Tag_>::value(sm, fine_index, coarse_index);
            std::cout << sm;
            SparseMatrixELL<DT1_> prolmat(sm);

            Prolongation<Tag_, methods::NONE>::value(defect_fine, defect_coarse, mask, prolmat);
            DenseVector<DT1_> ref_res(defect_fine.copy());
            Prolongation<Tag_, methods::PROLMAT>::value(defect_fine, defect_coarse, mask, prolmat);

            TEST_CHECK_EQUAL(defect_fine, ref_res);

        }
};
ProlongationPROLMATTest<tags::CPU, double> prolmat_cpu_double_test("CPU double", 81, 25);
