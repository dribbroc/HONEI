/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef BLASTEST_HH
#define BLASTEST_HH 1
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
extern "C" 
{
	#include "cblas.h"
}
#include "basetest.hh"
using namespace std;	
class BlasTest : public BaseTest
{
	public:		
		//BlasTest();
		BlasTest(string n);
		bool start_test();

};

#endif
