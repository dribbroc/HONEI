/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <string>
#include <exception>
#include <list>
#include <time.h>
#include <iostream>

using namespace std;
/**
 * Baseclass for all testingclasses
 */
class Benchmark
{
    protected:
        const std::string _id;
        double _start, _end;
		list<double> _benchlist;

    public:
        /**
		* Constructor
		* \param id the Benchmark
		*/
        Benchmark(const std::string & id);

        const std::string id() const;

        /// called by the benchmark framework to run the benchmark
        virtual void run() = 0;

		///generates a standard benchmark output
        void evaluate();

        ///generates a standard benchmark output with flop being the number of floating point operations of the benchmarked function
        void evaluate(int flop);
};

/**
 * Exception thrown by Benchmarks
 */
class BenchFailedException :
    public std::exception
{
    private:
        std::string _message;

    public:
        /**
         * Constructor.
         */
        BenchFailedException(const char * const function, const char * const file,
            const long line, const std::string & message) throw ()
        {
            _message = message;
        }

        /**
         * Destructor.
         */
        virtual ~BenchFailedException() throw ()
        {
        }

        /**
         * Description.
         */
        const char * what() const throw ()
        {
            return _message.c_str();
        }
};

/**
 * Benchmarks and adds result to list
 */
#define BENCHMARK(a) \
    do { \
        try { \
	    _start = clock(); \
	    a; \
	    _end = clock(); \
	    _benchlist.push_back(_end - _start); \
        } catch (const std::exception & test_e) { \
            throw BenchFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Benchmark threw exception."); \
        } catch (...) { \
            throw BenchFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Benchmark threw unknown exception."); \
        } \
    } while (false)

