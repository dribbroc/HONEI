/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/visual/engine_server.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>


using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class EngineServerTest :
    public BaseTest
{
    public:
        EngineServerTest(const std::string & type) :
            BaseTest("EngineServer test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            EngineServer<Tag_, DataType_> engine_server;
            engine_server.run();
            TEST_CHECK(true);
        }
};
//EngineServerTest<tags::CPU, float> engine_server_test_double("float");
