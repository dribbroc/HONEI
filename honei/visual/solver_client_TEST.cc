/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LibVisual C++ library. LibVisual is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibVisual is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/visual/solver_client.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>


using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class SolverClientTest :
    public BaseTest
{
    public:
        SolverClientTest(const std::string & type) :
            BaseTest("SoverClient test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<double> hf(4, 4);
            SolverClient<Tag_, DataType_> solver_client;

            solver_client.init("localhost", 4711, 12345);
            solver_client.do_step(hf);
            solver_client.quit_server();

            solver_client.init("localhost", 4711, 12345);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.restart_scenario();
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.quit_server();

            solver_client.init("localhost", 4711, 12345);
            solver_client.do_step(hf);
            solver_client.do_step(hf);
            solver_client.shutdown_server();
            TEST_CHECK(true);
        }
};
//SolverClientTest<tags::CPU, float> solver_client_test_double("float");
