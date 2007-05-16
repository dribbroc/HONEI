/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef TESTFACTORY_HH
#define TESTFACTORYT_HH 1
#include <string>
#include "basetest.hh"
#include "blastest.hh"
#include "liblatest.hh"
using namespace std;	
class TestFactory
{
	public:
	static BaseTest * new_test (const string & description);
};

#endif
