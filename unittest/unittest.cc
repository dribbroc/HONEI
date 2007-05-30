/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <unittest/unittest.hh>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <list>
#include <string>
#include <utility>

using namespace pg512;

class BaseTest;

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
    std::cout << "." << std::flush;
    if (! was_ok)
        throw TestFailedException(function, file, line, message);
}

TestFailedException::TestFailedException(const char * const function, const char * const file,
        const long line, const std::string & message) throw () :
    _message(pg512::stringify(file) + ":" + pg512::stringify(line) + ": in " +
            pg512::stringify(function) + ": " + message )
{
}

TestFailedException::~TestFailedException() throw ()
{
}

int main(int argc, char** argv) 
{
    bool quick = false;
  	if (argc==2 && stringify(argv[1])=="quick") quick=true;
    int result=EXIT_SUCCESS;

    for (TestList::Iterator i(TestList::instance()->begin_tests()),i_end(TestList::instance()->end_tests()) ; i != i_end ; ++i)
    {
        try
        {
            std::cout << (*i)->id() << ": " << std::flush;
            if (quick)
            {
                if ((*i)->is_quick_test()) (*i)->run();
            }
            else (*i)->run();
            std::cout << " PASSED" << std::endl;      
        }
        catch (TestFailedException & e)
        {
            std::cout << "FAILED: " << e.what() << std::endl;
            result = EXIT_FAILURE;
        }
    }

    return result;
}
