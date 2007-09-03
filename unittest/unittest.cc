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

#include <cstdlib>
#include <exception>
#include <iostream>
#include <list>
#include <string>
#include <utility>

using namespace honei;
using namespace tests;

class TestList
{
    private:
        static std::list<BaseTest *> _tests;

        TestList()
        {
        }

    public:
        typedef std::list<BaseTest*>::const_iterator Iterator;

        static TestList * instance()
        {
            static TestList result;

            return &result;
        }

        void register_test(BaseTest * const test)
        {
            _tests.push_back(test);
        }

        Iterator begin_tests() const
        {
            return _tests.begin();
        }

        Iterator end_tests() const
        {
            return _tests.end();
        }
};

std::list<BaseTest *> TestList::_tests;

BaseTest::BaseTest(const std::string & id) :
    _id(id)
{
    TestList::instance()->register_test(this);
}

const std::string BaseTest::id() const
{
    return _id;
}

void
BaseTest::check(const char * const function, const char * const file,
    const long line, bool was_ok, const std::string & message) const
{
    if (! was_ok)
        throw TestFailedException(function, file, line, message);
}

bool
BaseTest::is_quick_test() const
{
    return false;
}

QuickTest::QuickTest(const std::string & id) :
    BaseTest(id)
{
}

bool
QuickTest::is_quick_test() const
{
    return true;
}

TestFailedException::TestFailedException(const char * const function, const char * const file,
        const long line, const std::string & message) throw () :
    _message(honei::stringify(file) + ":" + honei::stringify(line) + ": in " +
            honei::stringify(function) + ": " + message )
{
}

TestFailedException::~TestFailedException() throw ()
{
}

int main(int argc, char** argv)
{
    int result(EXIT_SUCCESS);
    bool quick(false);

    if ((argc == 2) && (stringify(argv[1]) == "quick"))
            quick=true;

    for (TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
            i != i_end ; ++i)
    {
        CONTEXT("When running test case '" + (*i)->id() + "':");
        try
        {
            if (quick && (!(*i)->is_quick_test()) )
                continue;

            std::cout << (*i)->id() + ": \n";
            (*i)->run();
            std::cout << "PASSED \n";
        }
        catch (TestFailedException & e)
        {
            std::cout << "FAILED: " << std::endl << stringify(e.what()) << std::endl;
            result = EXIT_FAILURE;
        }
    }

    return result;
}
