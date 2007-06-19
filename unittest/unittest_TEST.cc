/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Andre Matuschek <andre@matuschek.org>
 * Copyright (c) 2007 David Gies      <david-gies@gmx.de> 
 * Copyright (c) 2007 Danny van Dyk   <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock   <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
 
#include <unittest/unittest.hh>
using namespace std;
using  namespace tests;

/*
 * Test class for the unittest framework itself
 */
class UnitTest : public BaseTest
{
    /// Constructor
	UnitTest(const std::string & id) : 
	    BaseTest(id)
	{
	}
	
	/// runs the tests
	virtual void run() const
	{
		TEST_CHECK(true);
		TEST_CHECK_EQUAL(1,1);
        TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
        TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
		TEST_CHECK_THROWS(string("0").at(10), exception);
	}
}unittest("UnitTest-test;


