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
#include <honei/lbm/solver_labswe.hh>
#include <unittest/unittest.hh>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class SolverLABSWETest :
    public BaseTest
{
    public:
        SolverLABSWETest(const std::string & type) :
            BaseTest("solver_labswe_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> dummy_m(10, 10);

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver
                (1.,1.,1., 1, 1, &dummy_m, &dummy_m, &dummy_m, &dummy_m);

            solver.do_preprocessing();
            solver.solve();
            TEST_CHECK(true);
        }

};
SolverLABSWETest<tags::CPU, float> solver_test_float("float");
