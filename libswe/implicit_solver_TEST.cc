/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
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

#include "implicit_solver.hh"
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <string>
#include <libswe/scenario.hh>

#include <sys/time.h>

using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;
using namespace boundaries;

template <typename Tag_, typename DataType_>
class ImplicitSolverCreationTest :
    public BaseTest
{
    public:
        ImplicitSolverCreationTest(const std::string & type) :
            BaseTest("implicit_solver_creation_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(2000);
            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);
            TEST_CHECK(true);
        }
};
ImplicitSolverCreationTest<tags::CPU, float> implicit_solver_creation_test_float("float");
