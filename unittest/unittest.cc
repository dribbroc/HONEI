/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifdef MOO
#include <unittest/unittest.hh>
using namespace std;

class BaseTest;

class TestList
{
    private:
        static list<BasteTest *> _tests;
        
        TestList()
        {
        }
        
    public:
        typedef libwrapiter::ForwardIterator<TestList, list<BaseTest*>::value_type> Iterator;

        static TestList * instance()
        {
            static TestList result;
            return &result;
        }
        
        void register_test(BaseTest * const test)
        {
            _tests.push_back(test);
        }
        
        Iterator begin_tests() const
        {
            return Iterator(_tests.begin());
        }
        
        Iterator end_tests() const
        {
            return Iterator(_tests.end());
        }
};

list<BaseTest *> TestList::_tests;

        BaseTest::BaseTest(const string & id) : _id(id)
        {
            TestList::instance()->register_test(this);
        }
        
        const string BaseTest::id() const
        {
            return _id;
        }

        void BaseTest::check(const char * const function, const char * const file,
 	        const long line, bool was_ok, const std::string & message) const
	 	{
	 	    std::cout << "." << std::flush;
	 	    if (! was_ok)
	 	        throw TestFailedException(function, file, line, message);
	 	}





int main(int argc, char** argv) 
{

    /*cout << "Hallo du!\n";
    cout << argc << " Argumente. Wert:  " << argv[1];
    cout << "\n\n";

    cout << "Wieviele Parameter ausgeben?" << endl;

    int anzahl;
    cin >> anzahl;
    cout << anzahl << endl;
    	
    for (int i = 0; (i < anzahl) && (i < argc); i++) {
    	cout << argv[i] << endl; 
    } */
    
    for (TestList::Iterator i(TestList::instance()->begin_tests()),i_end(TestList::instance()->end_tests()) ; i != i_end ; ++i)
    {
        try
        {
            (*i)->run();
        }
        catch (TestError & e)
        {
            
        }
    }
	return 0;
}
#endif





 
