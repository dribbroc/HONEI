/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifdef MOO
#include <unittest/unittest.hh>
using namespace std;

class UnitTest : public BaseTest
{
	UnitTest(const string & id) : BaseTest(id)
	{
	}
	
	virtual void run() const
	{
		//TEST_CHECK(true);
		//TEST_CHECK_EQUAL(1,1);
        //TEST_CHECK_STRINGIFY_EQUAL(4711, 4711);
        //TEST_CHECK_EQUAL_EPS(25,23,2.2);
		//TEST_CHECK_THROWS(string("0").at(10), exception);
	}
}unittest;
#endif

