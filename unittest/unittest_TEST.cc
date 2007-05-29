/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <unittest/unittest.hh>
using namespace std;

/*
 * Test class for the unittest framework itself
 */
class UnitTest : public BaseTest
{
    /// Constructor
	UnitTest(const std::string & id) : 
	    BaseTest(id)
	{
	}
	
	/// runs the tests
	virtual void run() const
	{
		TEST_CHECK(true);
		TEST_CHECK_EQUAL(1,1);
        TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
        TEST_CHECK_EQUAL_WITHIN_EPS(25,23,2.2);
		TEST_CHECK_THROWS(string("0").at(10), exception);
	}
}unittest("UnitTest-test;


