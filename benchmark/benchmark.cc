/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.hh>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <list>
#include <string>
#include <utility>

// using namespace pg512;
using namespace std;

class Benchmark;

class BenchmarkList
{
    private:
        static std::list<Benchmark *> _tests;

        BenchmarkList()
        {
        }

    public:
        typedef std::list<Benchmark*>::const_iterator Iterator;

        static BenchmarkList * instance()
        {
            static BenchmarkList result;

            return &result;
        }

        void register_test(Benchmark * const test)
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

std::list<Benchmark *> BenchmarkList::_tests;

Benchmark::Benchmark(const std::string & id) :
    _id(id)
{
    BenchmarkList::instance()->register_test(this);
}

void Benchmark::evaluate()
		{
			int x = 0;
			double a = 0;
			for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
			{
				a += *i;
				++x;
			}
			double b = (a/CLOCKS_PER_SEC);
			cout << "Runtime - total: " << b << "sec" << endl;
			_benchlist.sort();
			b = (_benchlist.front()/CLOCKS_PER_SEC);
			cout << "Runtime - lowest: " << b << "sec" << endl;
			b = (_benchlist.back()/CLOCKS_PER_SEC);
			cout << "Runtime - highest: " << b << "sec" << endl;
			b = ((a/x)/CLOCKS_PER_SEC);
			cout << "Runtime - average: " << b << "sec" << endl;
		}

const std::string Benchmark::id() const
{
    return _id;
}

int main(int argc, char** argv)
{
    int result=EXIT_SUCCESS;

    for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_tests()),i_end(BenchmarkList::instance()->end_tests()) ; i != i_end ; ++i)
    {
        try
        {
            std::cout << (*i)->id() << ": " << std::endl << std::endl;
            (*i)->run();
            std::cout << "Finished " << (*i)->id() << " succesfull!" << std::endl;
        }
        catch (TestFailedException & e)
        {
            std::cout << " FAILED" << std::endl;
			std::cout << e.what() << std::endl;
            result = EXIT_FAILURE;
        }
    }

    return result;
}