/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Andre Matuschek <andre@matuschek.org>
 * Copyright (c) 2007 David Gies      <david-gies@gmx.de>
 * Copyright (c) 2007 Danny van Dyk   <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Dirk Ribbrock   <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/unittest.hh>
#include <honei/util/instantiation_policy-impl.hh>
#ifdef HONEI_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <exception>
#include <iostream>
#include <list>
#include <string>
#include <utility>

using namespace honei;
using namespace tests;

class TestList :
        public InstantiationPolicy<TestList, Singleton>
{
    private:
        static std::list<BaseTest *> _tests;

    public:
        friend class InstantiationPolicy<TestList, Singleton>;
        typedef std::list<BaseTest*>::iterator Iterator;

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

        unsigned long size()
        {
            return _tests.size();
        }

        Iterator erase (Iterator i)
        {
            return _tests.erase(i);
        }

};

std::list<BaseTest *> TestList::_tests;

BaseTest::BaseTest(const std::string & id) :
    _id(id),
    _tag_name(tags::NONE::name)
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

void
BaseTest::register_tag(std::string tag_name)
{
    _tag_name = tag_name;
}

std::string
BaseTest::get_tag_name()
{
    return _tag_name;
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
#ifdef HONEI_MPI
    MPI_Init(&argc, &argv);
#endif

    int result(EXIT_SUCCESS);
    bool quick(false);
    bool sse(false);
    bool itanium(false);
    bool cuda(false);
    bool cell(false);
    bool mc(false);
    bool cpu(false);
    bool opencl(false);
    bool generic(false);
    bool all(true);
    if ((argc > 1) && (stringify(argv[1]) == "quick"))
    {
        quick = true;
    }

    if (argc > 2)
    {
        for (int index(2) ; index < argc ; ++index)
        {
            if (stringify(argv[index]) == "sse")
            {
                sse = true;
                all = false;
            }
            if (stringify(argv[index]) == "itanium")
            {
                itanium = true;
                all = false;
            }
            if (stringify(argv[index]) == "cuda")
            {
                cuda = true;
                all = false;
            }
            if (stringify(argv[index]) == "generic")
            {
                generic = true;
                all = false;
            }
            if (stringify(argv[index]) == "cell")
            {
                cell = true;
                all = false;
            }
            if (stringify(argv[index]) == "mc")
            {
                mc = true;
                all = false;
            }
            if (stringify(argv[index]) == "cpu")
            {
                cpu = true;
                all = false;
            }
            if (stringify(argv[index]) == "opencl")
            {
                opencl = true;
                all = false;
            }
        }
    }

    unsigned long list_size(0);

    for (TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
            i != i_end ; )
    {
            if (quick && (!(*i)->is_quick_test()) )
            {
                i = TestList::instance()->erase(i);
                continue;
            }
            if (!all)
            {
                if (((*i)->get_tag_name()=="sse") && !sse)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="generic") && !generic)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="itanium") && !itanium)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="cuda") && !cuda)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="mc-cuda") && (!mc && !cuda))
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="cell") && !cell)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="mc-sse") && (!mc && !sse))
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="mc-generic") && (!mc && !generic))
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="mc") && (!mc && !cpu))
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="cpu") && !cpu)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="opencl-cpu" || (*i)->get_tag_name()=="opencl-gpu") && !opencl)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
                if (((*i)->get_tag_name()=="none") && !all)
                {
                    i = TestList::instance()->erase(i);
                    continue;
                }
            }
            ++i;
            list_size++;
    }

    unsigned long iterator_index(1);
    for (TestList::Iterator i(TestList::instance()->begin_tests()), i_end(TestList::instance()->end_tests()) ;
            i != i_end ; )
    {
        CONTEXT("When running test case '" + (*i)->id() + "':");
        try
        {
            /*if (quick && (!(*i)->is_quick_test()) )
                continue;
            if (!all)
            {
                if (((*i)->get_tag_name()=="sse") && !sse)
                    continue;
                if (((*i)->get_tag_name()=="cuda") && !cuda)
                    continue;
                if (((*i)->get_tag_name()=="cell") && !cell)
                    continue;
                if (((*i)->get_tag_name()=="mc-sse") && (!mc && !sse))
                    continue;
                if (((*i)->get_tag_name()=="mc") && (!mc && !cpu))
                    continue;
                if (((*i)->get_tag_name()=="cpu") && !cpu)
                    continue;
                if (((*i)->get_tag_name()=="none") && !all)
                    continue;
            }*/

            std::cout << "(" << iterator_index << "/" << list_size << ") " << (*i)->id() + " [Backend: "
                << (*i)->get_tag_name() << "]" << std::endl;
            (*i)->run();
            std::cout << "PASSED" << std::endl;
        }
        catch (TestFailedException & e)
        {
            std::cout << "FAILED: " << std::endl << stringify(e.what()) << std::endl;
            result = EXIT_FAILURE;
        }
        catch (honei::Exception & e)
        {
            std::cout << "Caught exception:" << std::endl << e.message() << " " <<e.what() << std::endl;
            throw;
        }
        catch (std::exception & e)
        {
            std::cout << "Caught exception:" << std::endl << e.what() << std::endl;
            throw;
        }
        i = TestList::instance()->erase(i);
        iterator_index++;
    }

#ifdef HONEI_MPI
    MPI_Finalize();
#endif
    return result;
}
