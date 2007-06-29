/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <solver.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class RelaxSolverQuickTest :
    public QuickTest
{
    public:
        RelaxSolverQuickTest(const std::string & type) :
            QuickTest("relax_solver_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            //RelaxSolver<DataType_, DataType_, DataType_, DataType_, DataType_>  relax_solver;
                //= new RelaxSolver<DataType_, DataType_, DataType_, DataType_, DataType_>();
        }
};
RelaxSolverQuickTest<float> relax_solver_quick_test_float("float");
RelaxSolverQuickTest<double> relax_solver_quick_test_double("double");


