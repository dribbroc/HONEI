/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifdef MOO
#ifndef UNITTEST_HH
#define UNITTEST_HH

#include <iostream>
#include <string>
#include <exception>
#include <utility>
#include <list>
#include <libwraptiter/libwrapiter_forward_iterator.hh>
#include <libutil/stringify.hh>

using namespace std;

class BaseTest
{
    protected:
        const string _id;
    public:
        BaseTest(const string & id) : _id(id);
        
        const string id() const;
        
        
        virtual void run() const = 0;
        virtual void benchmark() const = 0;

        void check(const char * const function, const char * const file,
 	        const long line, bool was_ok, const std::string & message) const;	

    /**
 	     * Utility class used by TEST_CHECK_EQUAL.
 	     */
 	    struct TwoVarHolder
 	    {
 	        bool result;
 	        std::string s_a;
 	        std::string s_b;
 	
 	        template <typename T1_, typename T2_>
 	        TwoVarHolder(T1_ a, T2_ b) :
 	            result(a == b),
 	            s_a(stringify(a)),
 	            s_b(stringify(b))
 	        {
 	        }
 	    }; 	

        struct TestEps        
 	    {
 	        bool result;
 	        string s_a;
 	        string s_b;
 	        string s_diff;
 	
 	        template <typename T1_, typename T2_, typename T3_>
 	        TwoVarHolder(T1_ a, T2_ b, T3_ c) :
 	            result((a - b) < c),
 	            s_diff(stringify(a-b)),
 	            s_a(stringify(a)),
 	            s_b(stringify(b))
 	        {
 	        }
 	    }; 	
};

class TestFailedException : public exception
 	    {
 	        private:
 	            const string _message;
 	
 	        public:
 	            /**
 	             * Constructor.
 	             */
 	            TestFailedException(const char * const function, const char * const file,
 	                    const long line, const string & message) throw ();
 	
 	            /**
 	             * Destructor.
 	             */
 	            virtual ~TestFailedException() throw ();
 	
 	            /**
 	             * Description.
 	             */
 	            const char * what() const throw ()
 	            {
 	                return _message.c_str();
 	            }
 	    };

/**
 	 * Check that a == b.
 	 */
 	#define TEST_CHECK_EQUAL(a, b) \
 	    do { \
 	        try { \
 	            TwoVarHolder test_h(a, b); \
 	            check(__PRETTY_FUNCTION__, __FILE__, __LINE__, test_h.result, \
 	                    "Expected '" #a "' to equal '" + test_h.s_b + \
 	                    "' but got '" + test_h.s_a + "'"); \
 	        } catch (const TestFailedException &) { \
 	            throw; \
 	        } catch (const std::exception & test_e) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected exception " + \
 	                    " inside a TEST_CHECK_EQUAL block"); \
 	        } catch (...) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL block"); \
 	        } \
 	    } while (false)
 	
 	/**
 	 * Check that stringify(a) == stringify(b).
 	 */
 	#define TEST_CHECK_STRINGIFY_EQUAL(a, b) \
 	    do { \
          try { \
 	            std::string s_a(stringify(a)); \
 	            std::string s_b(stringify(b)); \
 	            check(__PRETTY_FUNCTION__, __FILE__, __LINE__, s_a == s_b, \
 	                    "Expected '" #a "' to equal '" + s_b + \
 	                    "' but got '" + s_a + "'"); \
 	        } catch (const TestFailedException &) { \
 	            throw; \
 	        } catch (const exception & test_e) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected exception " + \
 	                    " inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
 	        } catch (...) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected unknown exception inside a TEST_CHECK_STRINGIFY_EQUAL block"); \
 	        } \
 	    } while (false)
 	
 	/**
 	 * Check that a is true.
 	 */
 	#define TEST_CHECK(a) \
 	    do { \
 	        try { \
 	            check(__PRETTY_FUNCTION__, __FILE__, __LINE__, a, \
 	                    "Check '" #a "' failed"); \
 	        } catch (const TestFailedException &) { \
 	            throw; \
 	        } catch (const std::exception & test_e) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected exception "  + \
 	                    " inside a TEST_CHECK block"); \
 	        } catch (...) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected unknown exception inside a TEST_CHECK block"); \
 	        } \
 	    } while (false)
 	
 	/**
 	 * Check that a throws an exception of type b.
 	 */
 	#define TEST_CHECK_THROWS(a, b) \
 	    do { \
 	        try { \
 	            try { \
	                a; \
 	                check(__PRETTY_FUNCTION__, __FILE__, __LINE__, false, \
 	                        "Expected exception of type '" #b "' not thrown"); \
 	            } catch (b &) { \
 	                TEST_CHECK(true); \
 	            } \
 	        } catch (const TestFailedException &) { \
 	            throw; \
 	        } catch (const std::exception & test_e) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected exception "  + \
 	                    " inside a TEST_CHECK_THROWS block"); \
 	        } catch (...) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected unknown exception inside a TEST_CHECK_THROWS block"); \
 	        } \
 	    } while (false)

/**
 	 * Check that a -b < epsilon.
 	 */
 	#define TEST_CHECK_EQUAL_EPS(a, b, eps) \
 	    do { \
 	        try { \
 	            TestEps test_eps(a, b, eps); \
 	            check(__PRETTY_FUNCTION__, __FILE__, __LINE__, test_eps.result, \
 	                    "Expected '" #a "' - '" + #b + \
 	                    "' < '" + eps + "'but was'"+ test_eps.s_diff); \
 	        } catch (const TestFailedException &) { \
 	            throw; \
 	        } catch (const std::exception & test_e) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected exception " + \
 	                    " inside a TEST_CHECK_EQUAL block"); \
 	        } catch (...) { \
 	            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
 	                    "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL block"); \
 	        } \
 	    } while (false)


#endif
