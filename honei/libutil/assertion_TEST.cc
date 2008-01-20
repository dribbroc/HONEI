/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/assertion.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using namespace tests;

class AssertionTest :
    public QuickTest
{
    public:
        AssertionTest() :
            QuickTest("assertion_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_THROWS(ASSERT(false, "Should throw!"), Assertion);

            bool no_exception_thrown(true);
            try
            {
                ASSERT(true, "Shouldn't throw!");
            }
            catch (...)
            {
                no_exception_thrown = false;
            }
            TEST_CHECK(no_exception_thrown);
        }
} assertion_test;
