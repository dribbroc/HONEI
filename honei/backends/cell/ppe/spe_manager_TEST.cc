/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/util/mutex.hh>
#include <honei/util/lock.hh>
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
            register_tag(tags::Cell::name);
        }

        virtual void run() const
        {
            try
            {
                /*unsigned long count(0);
                for (SPEManager::Iterator i(SPEManager::instance()->begin()),
                        i_end(SPEManager::instance()->end()) ; i != i_end ; ++i, ++count)
                {
                    std::cout << "SPE id is '" << i->id() << "'" << std::endl;
                }
                std::cout << "Number of SPEs is " << count << std::endl;*/
                std::cout << "Number of SPEs is " << SPEManager::instance()->spe_count() << std::endl;
                TEST_CHECK(true);
            }
            catch (...)
            {
                TEST_CHECK(false);
            }
        }
} spe_manager_iteration_test;
