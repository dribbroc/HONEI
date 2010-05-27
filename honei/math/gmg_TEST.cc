/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.tu-dortmund.de>
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

#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <cmath>
#include <honei/math/gmg.hh>
#include <honei/math/methods.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class GMGTest:
    public BaseTest
{
    public:
        GMGTest(const std::string & tag) :
            BaseTest("GMG test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            GMGInfo<DT1_> info;

            GMG<Tag_, methods::NONE>::value(info);
        }
};
GMGTest<tags::CPU, float> cpu_gmg_test_float("float");

