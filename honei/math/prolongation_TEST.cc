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
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/la/dense_matrix.hh>
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
            unsigned long N_fine = 1089;
            unsigned long width_fine = (unsigned long)sqrt((double)N_fine);
            unsigned long N_coarse = 289;
            unsigned long width_coarse = (unsigned long)sqrt((double)N_coarse);

            DenseVector<DT1_> fine(N_fine, DT1_(4711));
            DenseVector<DT1_> fine_ref(N_fine);
            DenseVector<DT1_> coarse(N_coarse);
            for (unsigned long i(0) ; i < coarse.size() ; ++i)
            {
                coarse[i] = DT1_(i % 1000);
            }
            DenseVector<unsigned long> mask(8);
            for(unsigned long i(0) ; i < 8 ; ++i)
            {
                //if(i % 2  == 0)
                    mask[i] = 2;
            }

            Prolongation<Tag_>::value(fine, coarse, mask);
            Prolongation<tags::CPU>::value(fine_ref, coarse, mask);
            fine.lock(lm_read_only);
            fine_ref.lock(lm_read_only);
            TEST_CHECK_EQUAL(fine, fine_ref);

            /*DenseMatrix<DT1_> grid_fine(width_fine, width_fine, DT1_(0));
            for(unsigned long i(0) ; i < width_fine ; ++i)
            {
                for(unsigned long j(0); j < width_fine; ++j)
                {
                    grid_fine(i,j) = fine[i * j + j];
                }
            }
            std::cout << grid_fine << std::endl;*/
            fine.unlock(lm_read_only);
            fine_ref.unlock(lm_read_only);

            TEST_CHECK(true);
        }
};
ProlongationTest<tags::CPU, float> prolongate_test("float");
#ifdef HONEI_CUDA
ProlongationTest<tags::GPU::CUDA, float> cuda_prolongate_test("float");
#endif
