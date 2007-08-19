/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <libutil/mutex.hh>
#include <libutil/lock.hh>
#include <libutil/spe_manager.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace pg512;
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

namespace
{
    class TestTask
    {
        private:
            static DeviceId v;
            Mutex * const m;

        public:
            TestTask() :
                m(new Mutex)
            {
            }

            virtual void operator() (const SPE & spe)
            {
                DeviceId id(spe.id());
                Lock l(*m);

                v += id;
            }

            inline DeviceId value() const
            {
                Lock l(*m);

                return v;
            }
    };

    DeviceId TestTask::v = 0;
}

class SPEManagerThreadingTest :
    public QuickTest
{
    public:
        SPEManagerThreadingTest() :
            QuickTest("spe_manager_threading_test")
        {
        }

        virtual void run() const
        {
            TestTask t;
            SPETask st(t);
            unsigned count(std::distance(SPEManager::instance()->begin(), SPEManager::instance()->end())),
                     repetitions(5);
            bool unequal;

            SPEManager::instance()->dispatch(st);

            while ((unequal = ((count * (count - 1) / 2) != t.value())) && (repetitions > 0))
            {
                sleep(1);
                --repetitions;
            }

            TEST_CHECK(! unequal);
        }
} spe_manager_threading_test;
