/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef UNITTEST_GUARD_UNITTEST_HH
#define UNITTEST_GUARD_UNITTEST_HH 1

#include <libutil/stringify.hh>

#include <string>
#include <exception>

/**
 * Baseclass for all testingclasses
 */
class BaseTest
{
    protected:
        const std::string _id;

    public:
        /**
         * Constructor
         * \param id the test-name
         */
        BaseTest(const std::string & id);

        const std::string id() const;

        /// Utility method used bei TEST_CHECK_*
        void check(const char * const function, const char * const file,
                const long line, bool was_ok, const std::string & message) const;

        /// called by the unittest framework to run the tests
        virtual void run() const = 0;
        
        /// returns true if the test is a "quick test"
        virtual bool is_quick_test()
        {
            return false;
        }
        
        /// returns testclass id string
        virtual std::string get_id()
        {
            return _id;
        }

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

        /**
         * Utility class used by TEST_CHECK_EQUAL_WITHIN_EPS.
         */
        struct WithinEpsCalculator
        {
            bool result;
            std::string s_a;
            std::string s_b;
            std::string s_diff;

            template <typename T1_, typename T2_, typename T3_>
            WithinEpsCalculator(T1_ a, T2_ b, T3_ c) :
                s_a(stringify(a)),
                s_b(stringify(b))
            {
				if (a>=b) 
				{
					result = ((a - b) <= c);
					s_diff = stringify(a-b);
				}
				else 
				{
					result = ((b - a) <= c);
					s_diff = stringify(b-a);
				}
            }
        };
};

/**
 * Baseclass for all testingclasses, marked as quick tests
 */
class QuickTest : public BaseTest
{
    public:
        QuickTest(const std::string & id): BaseTest(id) {}
        
        virtual bool is_quick_test()
        {
            return true;
        }
};

/**
 * Exception thrown by the check method in BaseTest
 */
class TestFailedException :
    public std::exception
{
    private:
        const std::string _message;

    public:
        /**
         * Constructor.
         */
        TestFailedException(const char * const function, const char * const file,
            const long line, const std::string & message) throw ();


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
                    _id + "\n" +  "Expected '" #a "' to equal '" + test_h.s_b + \
                    "' but got '" + test_h.s_a + "'"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ pg512::stringify(test_e.what()) + \
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
                    _id + "\n" +  "Expected '" #a "' to equal '" + s_b + \
                    "' but got '" + s_a + "'"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception  "+ pg512::stringify(test_e.what()) + \
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
                    _id + "\n" +  "Check '" #a "' failed"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ pg512::stringify(test_e.what()) + \
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
                        _id + "\n" +  "Expected exception of type '" #b "' not thrown"); \
            } catch (b &) { \
                TEST_CHECK(true); \
            } \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ pg512::stringify(test_e.what()) + \
                    " inside a TEST_CHECK_THROWS block"); \
        } catch (...) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected unknown exception inside a TEST_CHECK_THROWS block"); \
        } \
    } while (false)

/**
 * Check that a - b < epsilon.
 */
#define TEST_CHECK_EQUAL_WITHIN_EPS(a, b, eps) \
    do { \
        try { \
            WithinEpsCalculator calc(a, b, eps); \
            check(__PRETTY_FUNCTION__, __FILE__, __LINE__, calc.result,  \
                _id + "\n" + "Expected '|" #a " - " #b \
                "|' < '" + stringify(eps) + "' but was '" + calc.s_diff +"'"); \
        } catch (const TestFailedException & test_e) { \
            throw;  \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception  "+ pg512::stringify(test_e.what()) + \
                    " inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
        } catch (...) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
        } \
    } while (false)

#endif
