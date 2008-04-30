/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/lbm/tags.hh>
#include <unittest/unittest.hh>
#include <honei/lbm/equilibrium_distribution.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class EqDisLABSWETest :
    public BaseTest
{
    public:
        EqDisLABSWETest(const std::string & type) :
            BaseTest("eq_dis_labswe_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> h(1000ul, 1000ul, DataType_(1.23456));
            DenseMatrix<DataType_> u(1000ul, 1000ul, DataType_(1.23456));
            DenseMatrix<DataType_> v(1000ul, 1000ul, DataType_(1.23456));
            DataType_ g(9.81);
            DataType_ e(1.);
            DataType_ e_u(2.);
            DataType_ e_v(2.);

            DenseMatrix<DataType_> result_1(1000ul, 1000ul);
            DenseMatrix<DataType_> result_2(1000ul, 1000ul);
            DenseMatrix<DataType_> result_3(1000ul, 1000ul);

            EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_0>::
                value(result_1, h, g , e);

            EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_ODD>::
                value(result_2, h, u, v, g, e, e_u, e_v);

            EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_EVEN>::
                value(result_3, h, u, v, g, e, e_u, e_v);

            for(unsigned long i(0); i < 1000; ++i)
            {
                for(unsigned long j(0); j < 1000; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(result_1(i,j), (1.23456 - ((5. * 9.81 * 1.23456 * 1.23456) / (6.)) - ((2. * 1.23456) / (3.))), std::numeric_limits<DataType_>::epsilon() * 7.);
                    TEST_CHECK_EQUAL_WITHIN_EPS(result_2(i,j), ((9.81 * 1.23456 * 1.23456) / 6. + ((1.23456 / 3.) * 2. * 1.23456) + ((1.23456 / 2.) * 2. * 1.23456 * 1.23456) - ((1.23456 / 6.) * 1.23456 * 1.23456)), std::numeric_limits<DataType_>::epsilon() * 7.);
                    TEST_CHECK_EQUAL_WITHIN_EPS(result_3(i,j), ((9.81 * 1.23456 * 1.23456) / 24. + ((1.23456 / 12.) * 2. * 1.23456) + ((1.23456 / 8.) * 2. * 1.23456 * 1.23456) - ((1.23456 / 24.) * 1.23456 * 1.23456)), std::numeric_limits<DataType_>::epsilon() * 7.);
                }
            }
        }
};
EqDisLABSWETest<tags::CPU, double> source_test_double("CPU double");
EqDisLABSWETest<tags::CPU, float> source_test_float("CPU float");
EqDisLABSWETest<tags::CPU::MultiCore, double> source_test_double_mc("MC double");
EqDisLABSWETest<tags::CPU::MultiCore, float> source_test_float_mc("MC float");
#ifdef HONEI_SSE
EqDisLABSWETest<tags::CPU::SSE, double> source_test_double_sse("SSE double");
EqDisLABSWETest<tags::CPU::SSE, float> source_test_float_sse("SSE float");
EqDisLABSWETest<tags::CPU::MultiCore::SSE, double> source_test_double_mc_sse("MCSSE double");
EqDisLABSWETest<tags::CPU::MultiCore::SSE, float> source_test_float_mc_sse("MCSSE float");
#endif
