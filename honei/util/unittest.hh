/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Andre Matuschek <andre@matuschek.org>
 * Copyright (c) 2007 David Gies <david-gies@gmx.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef UNITTEST_GUARD_UNITTEST_HH
#define UNITTEST_GUARD_UNITTEST_HH 1

#include <honei/util/stringify.hh>
#include <honei/util/exception.hh>
#include <honei/util/tags.hh>

#include <string>
#include <exception>

/**
 * \file
 *
 * Implementation of BaseTest and related classes.
 *
 * \ingroup tests
 **/
using namespace honei;

namespace tests
{

    /**
     * Baseclass for all testingclasses
     * \ingroup tests
     */
    class BaseTest
    {
        protected:
            const std::string _id;
            std::string _tag_name;

        public:
            /**
             * Constructor.
             *
             * \param id The testcase's id string.
             */
            BaseTest(const std::string & id);

            /// Destructor.
            virtual ~BaseTest() {}

            /// Returns our id string.
            virtual const std::string id() const;

            /// Utility method used bei TEST_CHECK_*
            virtual void check(const char * const function, const char * const file,
                    const long line, bool was_ok, const std::string & message) const;

            /**
             * Runs the test case.
             *
             * Called by unittest framework only.
             * \ingroup tests
             */
            virtual void run() const = 0;

            /// Returns whether we are a quick-test.
            virtual bool is_quick_test() const;

            /// Register our target platform.
            virtual void register_tag(std::string tag_name);

            /// Returns our target platform.
            virtual std::string get_tag_name();

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
                    s_a(),
                    s_b()
                {
                    if (!result)
                    {
                        s_a = stringify(a);
                        s_b = stringify(b);
                    }
                }
            };

            /**
             * Utility class used by TEST_CHECK_NOT_EQUAL.
             */
            struct TwoVarHolder2
            {
                bool result;
                std::string s_a;
                std::string s_b;

                template <typename T1_, typename T2_>
                TwoVarHolder2(T1_ a, T2_ b) :
                    result(a == b),
                    s_a(),
                    s_b()
                {
                    if (result)
                    {
                        s_a = stringify(a);
                        s_b = stringify(b);
                    }
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
                    s_a(),
                    s_b()
                {
                    if (a >= b)
                    {
                        result = ((a - b) <= c);
                        if (!result)
                        {
                            s_diff = stringify(a - b);
                            s_a = stringify(a);
                            s_b = stringify(b);
                        }
                    }
                    else
                    {
                        result = ((b - a) <= c);
                        if (!result)
                        {
                            s_diff = stringify(b - a);
                            s_a = stringify(a);
                            s_b = stringify(b);
                        }
                    }
                }
            };
    };

    /**
     * Abstract Baseclass for all tagged test classes.
     * \ingroup tests
     */
    template <typename Tag_>
    class TaggedTest : public BaseTest
    {
        public:
            /**
             * Constructor.
             *
             * \param id The testcase's id string.
             */
            TaggedTest(const std::string & id):
                BaseTest(id)
            {
                _tag_name = Tag_::name;
            };

            /// Destructor.
            virtual ~TaggedTest() {}

            /// Returns whether we are a quick-test.
            virtual bool is_quick_test() const
            {
                return false;
            };
    };

    /**
     * Abstract Baseclass for all tagged quick test classes.
     * \ingroup tests
     */
    template <typename Tag_>
    class QuickTaggedTest : public BaseTest
    {
        public:
            /**
             * Constructor.
             *
             * \param id The testcase's id string.
             */
            QuickTaggedTest(const std::string & id):
                BaseTest(id)
            {
                _tag_name = Tag_::name;
            }

            /// Destructor.
            virtual ~QuickTaggedTest() {}

            /// Returns whether we are a quick-test.
            virtual bool is_quick_test() const
            {
                return true;
            }
    };

    /**
     * Abstract Baseclass for all untagged quick test classes.
     * \ingroup tests
     */
    class QuickTest : public BaseTest
    {
        public:
            /**
             * Constructor.
             *
             * \param id The testcase's id string.
             */
            QuickTest(const std::string & id):
                BaseTest(id)
            {
            }

            /// Destructor.
            virtual ~QuickTest() {}

            /// Returns whether we are a quick-test.
            virtual bool is_quick_test() const
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
}
/**
 * Check that a == b.
 */
#define TEST_CHECK_EQUAL(a, b) \
    do { \
        try { \
            BaseTest::TwoVarHolder test_h(a, b); \
            this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, test_h.result, \
                    this->_id + "\n" +  "Expected '" #a "' to equal \n'" + test_h.s_b + \
                    "'\nbut got\n'" + test_h.s_a + "'"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ honei::stringify(test_e.what()) + \
                    " inside a TEST_CHECK_EQUAL block"); \
        } catch (...) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL block"); \
        } \
    } while (false)

/**
 * Check that a != b.
 */
#define TEST_CHECK_NOT_EQUAL(a, b) \
    do { \
        try { \
            BaseTest::TwoVarHolder2 test_h(a, b); \
            this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, !test_h.result, \
                    this->_id + "\n" +  "Expected '" #a "' that is'" + test_h.s_a + \
                    "' to equal not '" + test_h.s_b + "'"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ honei::stringify(test_e.what()) + \
                    " inside a TEST_CHECK_NOT_EQUAL block"); \
        } catch (...) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected unknown exception inside a TEST_CHECK_NOT_EQUAL block"); \
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
            this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, s_a == s_b, \
                    this->_id + "\n" +  "Expected '" #a "' to equal '" + s_b + \
                    "'\nbut got\n'" + s_a + "'"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception  "+ honei::stringify(test_e.what()) + \
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
            this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, a, \
                    this->_id + "\n" +  "Check '" #a "' failed"); \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ honei::stringify(test_e.what()) + \
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
                this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, false, \
                        this->_id + "\n" +  "Expected exception of type '" #b "' not thrown"); \
            } catch (b &) { \
                TEST_CHECK(true); \
            } \
        } catch (const TestFailedException &) { \
            throw; \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception "+ honei::stringify(test_e.what()) + \
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
            BaseTest::WithinEpsCalculator calc(a, b, eps); \
            this->check(__PRETTY_FUNCTION__, __FILE__, __LINE__, calc.result,  \
                this->_id + "\n" + "Expected '|" #a "(" + stringify(a) + ") - " #b \
                "(" + stringify(b) + ")|' < '" + stringify(eps) + "' but was '" + calc.s_diff +"'"); \
        } catch (const TestFailedException & test_e) { \
            throw;  \
        } catch (const std::exception & test_e) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected exception  "+ honei::stringify(test_e.what()) + \
                    " inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
        } catch (...) { \
            throw TestFailedException(__PRETTY_FUNCTION__, __FILE__, __LINE__, \
                    "Test threw unexpected unknown exception inside a TEST_CHECK_EQUAL_WITHIN_EPS block"); \
        } \
    } while (false)

/**
 * Run the given test with pre- and postprocessing
 */
#define TEST(pre, test, post) \
    do { \
        pre; \
        try { \
            test; \
        } catch (const TestFailedException & test_e) { \
            post; \
            throw; } \
        post; \
    } while (false)

#endif
