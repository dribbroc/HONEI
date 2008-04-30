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
#include <honei/lbm/collide_stream.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class CollideStreamLABSWETest :
    public BaseTest
{
    public:
        CollideStreamLABSWETest(const std::string & type) :
            BaseTest("collideandstream_labswe_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> test(1000ul, 1000ul, DataType_(1.23456));
            DataType_ e_x(1.);
            DataType_ e_y(1.);
            DataType_ tau(1.);

            DenseMatrix<DataType_> result(1000ul, 1000ul);

            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_0>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_1>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_2>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_3>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_4>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_5>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_6>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_7>::
                value(result, test, test, test, test, e_x, e_y, tau);
            CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_PERIODIC, lbm_lattice_types::D2Q9::DIR_8>::
                value(result, test, test, test, test, e_x, e_y, tau);
            TEST_CHECK(true);
        }
};
CollideStreamLABSWETest<tags::CPU, float> source_test_float("CPU float");
CollideStreamLABSWETest<tags::CPU, double> source_test_double("CPU double");
CollideStreamLABSWETest<tags::CPU::MultiCore, float> source_test_float_mc("MC float");
CollideStreamLABSWETest<tags::CPU::MultiCore, double> source_test_double_mc("MC double");
#ifdef HONEI_SSE
CollideStreamLABSWETest<tags::CPU::SSE, float> source_test_float_sse("SSE float");
CollideStreamLABSWETest<tags::CPU::SSE, double> source_test_double_sse("SSE double");
CollideStreamLABSWETest<tags::CPU::MultiCore::SSE, float> source_test_float_mc_sse("MCSSE float");
CollideStreamLABSWETest<tags::CPU::MultiCore::SSE, double> source_test_double_mc_sse("MCSSE double");
#endif
