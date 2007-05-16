/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef LIBLATEST_HH
#define LIBLATEST_HH 1
#include <iostream>
#include "basetest.hh"
using namespace std;	
class LiblaTest : public BaseTest
{
	public:		
		//BlasTest();
		LiblaTest(string n);
		bool start_test();

};

#endif
