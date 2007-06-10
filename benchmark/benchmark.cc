/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "benchmark.hh"

#include <cstdlib>
#include <utility>
#include <fstream>

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
            double min = _benchlist.front();
            double max = _benchlist.front();
			for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
			{
				a += *i;
				++x;
                if (*i < min) min = *i;
                else if (*i > max) max = *i;
			}
			cout << "Function Calls: " << x << endl;
			double b = (a/CLOCKS_PER_SEC);
			cout << "Runtime - total: " << b << "sec" << endl;
			double c = (min/CLOCKS_PER_SEC);
			cout << "Runtime - lowest: " << c << "sec" << endl;
			double d = (max/CLOCKS_PER_SEC);
			cout << "Runtime - highest: " << d << "sec" << endl;
			double e = ((a/x)/CLOCKS_PER_SEC);
			cout << "Runtime - average: " << e << "sec" << endl;
			ofstream ofs("BenchmarkOut.txt", ios_base::out | ios_base::app);
			if (!ofs)
                cout << "Can't write to file!" << endl;
			else
			{
                ofs << _id  << " - " << __DATE__ << " - " << __TIME__ << endl << endl;
                ofs << "Result:"<< endl;
                ofs << "Function Calls: " << x << endl;
                ofs << "Runtime - total: " << b << "sec" << endl;
                ofs << "Runtime - lowest: " << c << "sec" << endl;
                ofs << "Runtime - highest: " << d << "sec" << endl;
                ofs << "Runtime - average: " << e << "sec" << endl;
                ofs << endl << endl << endl;
			}

		}

void Benchmark::evaluate(int flop)
		{
			int x = 0;
			double a = 0;
            double min = _benchlist.front();
            double max = _benchlist.front();
			for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
			{
				a += *i;
				++x;
                if (*i < min) min = *i;
                else if (*i > max) max = *i;
			}
			cout << "Function Calls: " << x << endl;
			double b = (a/CLOCKS_PER_SEC);
			cout << "Runtime - total: " << b << "sec" << endl;
			double c = (min/CLOCKS_PER_SEC);
			cout << "Runtime - lowest: " << c << "sec" << endl;
			double d = (max/CLOCKS_PER_SEC);
			cout << "Runtime - highest: " << d << "sec" << endl;
			double e = ((a/x)/CLOCKS_PER_SEC);
			cout << "Runtime - average: " << e << "sec" << endl;
			double f = 0;
			if (b > 0)
                f = ((x/b)*flop);
            if (f > 100000000)
                cout << (f/1000000000) << " GFLOPS" << endl;
            else
                if (f > 100000)
                    cout << (f/1000000) << " MFLOPS" << endl;
                else
                    if (f > 100)
                        cout << (f/1000) << " KFLOPS" << endl;
                    else
                        cout << (f) << " FLOPS" << endl;
			ofstream ofs("BenchmarkOut.txt", ios_base::out | ios_base::app);
			if (!ofs)
                cout << "Can't write to file!" << endl;
			else
			{
                ofs << _id  << " - " << __DATE__ << " - " << __TIME__ << endl << endl;
                ofs << "Result:"<< endl;
                ofs << "Function Calls: " << x << endl;
                ofs << "Runtime - total: " << b << "sec" << endl;
                ofs << "Runtime - lowest: " << c << "sec" << endl;
                ofs << "Runtime - highest: " << d << "sec" << endl;
                ofs << "Runtime - average: " << e << "sec" << endl;
                if (f > 100000000)
                    ofs << (f/1000000000) << " GFLOPS" << endl;
                else
                    if (f > 100000)
                        ofs << (f/1000000) << " MFLOPS" << endl;
                    else
                        if (f > 100)
                            ofs << (f/1000) << " KFLOPS" << endl;
                        else
                            ofs << (f) << " FLOPS" << endl;
                ofs << endl << endl << endl;
			}
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
        catch (BenchFailedException & e)
        {
            std::cout << " FAILED" << std::endl;
			std::cout << e.what() << std::endl;
            result = EXIT_FAILURE;
        }
    }

    return result;
}
