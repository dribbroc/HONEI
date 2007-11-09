/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <string> 
#include <exception>
#include <list>
#include <time.h>
#include <sys/time.h>
#include <iostream>

struct BenchmarkInfo
{
    unsigned long flops;
    unsigned long load;
    unsigned long store;
    
    BenchmarkInfo operator+(const BenchmarkInfo a) 
    {
        flops += a.flops;
        load += a.load;
        store += a.store;
        return *this;
    }

    BenchmarkInfo operator*(int a) 
    {
        flops *= a;
        load *= a;
        store *= a;
        return *this;
    }
};

using namespace std;
/**
 * Baseclass for all testingclasses
 */
class Benchmark
{
    protected:
        const std::string _id;
        timeval _start, _end;
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

        ///generates a benchmark output with FLOPS and data transfer rate information.
        void evaluate(BenchmarkInfo info);
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
        gettimeofday(&_start, 0); \
        a; \
        gettimeofday(&_end, 0); \
        _benchlist.push_back(((double)_end.tv_sec + ((double)_end.tv_usec/1000000)) - ((double)_start.tv_sec + ((double)_start.tv_usec/1000000))); \
        } catch (const std::exception & test_e) { \
            throw BenchFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Benchmark threw exception."); \
        } catch (...) { \
            throw BenchFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Benchmark threw unknown exception."); \
        } \
    } while (false)

