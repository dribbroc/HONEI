/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef TESTDISPATCHERT_HH
#define TESTDISPATCHER_HH 1
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include "basetest.hh"
#include "blastest.hh"
#include "liblatest.hh"
using namespace std;
/**
 *  A TestDispatcher dispatches a specified test to various subclasses of BaseTest
 */
class TestDispatcher
{
	public:
    /// runs the specified test and returns if it has failed or succeded
	bool start_test(int typ ,BaseTest * instances[], int instances_size);
    
    private:
    /// generate a random dense vector with given size
    //Vector * _random_dense_vector(unsigned long size);
    /// generate a random dense matrix with given size
    //Matrix * _random_dense_matrix(unsigned long columns, unsigned long rows);
    /// perform a read-write test on a vector
    //bool _vector_rw(BaseTest * instances[],int instances_size, Vector vector , int operation_count);
    /// perform a read-write test on a matrix
    //bool _matrix_rw(BaseTest * instances[],int instances_size,Matrix matrix, int operation_count);


};
#endif
