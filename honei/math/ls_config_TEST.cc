/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/unittest.hh>
#include <honei/math/ls_config.hh>

using namespace honei;
using namespace tests;

class LSConfigQuickTest :
    public QuickTest
{
    private:
        std::string _file;

    public:
        LSConfigQuickTest(std::string file) :
            QuickTest("ls_config_test: " + file),
            _file(file)
        {
        }

        virtual void run() const
        {
            ls_config::LSConfig::parse_file(_file);
        }
};
LSConfigQuickTest ls_config_test("lssolver.conf");
LSConfigQuickTest ls_config_test2("lssolver_simple.conf");
