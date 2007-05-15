/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef BLASTEST_HH
#define BLASTEST_HH 1
#include "basetest.hh"
using namespace std;	
class BlasTest : BaseTest
{
	public:		
		BlasTest();
		BlasTest(string n);
		bool starteTest();

};

#endif
