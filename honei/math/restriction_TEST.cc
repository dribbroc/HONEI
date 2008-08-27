/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/restriction.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/la/dense_matrix.hh>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class RestrictionTest:
    public BaseTest
{
    public:
        RestrictionTest(const std::string & tag) :
            BaseTest("Restriction test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long N_fine = 289;
            unsigned long width_fine = (unsigned long)sqrt((double)N_fine);
            unsigned long N_coarse = 81;
            unsigned long width_coarse = (unsigned long)sqrt((double)N_coarse);

            DenseVector<DT1_> fine(N_fine, DT1_(2));
            DenseVector<DT1_> coarse(N_coarse, DT1_(1));
            coarse[N_fine/2] = DT1_(4);
            DenseVector<unsigned long> mask(8);
            for(unsigned long i(0) ; i < 8 ; ++i)
            {
                if(i % 2  == 0)
                    mask[i] = 2;
            }

            Restriction<Tag_>::value(coarse, fine, mask);
            coarse.lock(lm_read_only);

            DenseMatrix<DT1_> grid_coarse(width_coarse, width_coarse, DT1_(0));
            for(unsigned long i(0) ; i < width_coarse ; ++i)
            {
                for(unsigned long j(0); j < width_coarse; ++j)
                {
                    grid_coarse(i,j) = coarse[i*j + j];
                }
            }
            coarse.unlock(lm_read_only);

            std::cout << grid_coarse << std::endl;
            TEST_CHECK(true);
        }
};
RestrictionTest<tags::CPU, float> restriction_test("float");
RestrictionTest<tags::GPU::CUDA, float> cuda_restriction_test("float");
