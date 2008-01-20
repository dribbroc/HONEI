/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/mutex.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/spe_manager.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace honei;
using namespace tests;

class SPEManagerIterationTest :
    public QuickTest
{
    public:
        SPEManagerIterationTest() :
            QuickTest("spe_manager_iteration_test")
        {
        }

        virtual void run() const
        {
            try
            {
                unsigned long count(0);
                for (SPEManager::Iterator i(SPEManager::instance()->begin()),
                        i_end(SPEManager::instance()->end()) ; i != i_end ; ++i, ++count)
                {
                    std::cout << "SPE id is '" << i->id() << "'" << std::endl;
                }
                std::cout << "Number of SPEs is " << count << std::endl;
                TEST_CHECK(true);
            }
            catch (...)
            {
                TEST_CHECK(false);
            }
        }
} spe_manager_iteration_test;
