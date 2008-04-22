#ifndef LIBUTIL_GUARD_BENCHMARK_INFO_HH
#define LIBUTIL_GUARD_BENCHMARK_INFO_HH 1

#include <list>
#include <iostream>

struct BenchmarkInfo
{
    unsigned long long flops;
    unsigned long long load;
    unsigned long long store;
    std::list<unsigned long> size;
    double scale;
    std::string scaleinfo;

    BenchmarkInfo() :
        flops(0),
        load(0),
        store(0),
        scale(1),
        scaleinfo(" ")
    {
    }
    
    BenchmarkInfo operator+(const BenchmarkInfo a) 
    {
        flops += a.flops;
        load += a.load;
        store += a.store;
        return *this;
    }

    BenchmarkInfo operator*(unsigned long a) 
    {
        flops *= a;
        load *= a;
        store *= a;
        return *this;
    }
};

#endif
