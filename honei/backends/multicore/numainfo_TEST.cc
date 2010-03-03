/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/backends/multicore/numainfo.hh>
#include <unittest/unittest.hh>

using namespace honei::mc;
using namespace tests;

class NumaInfoTest :
    public BaseTest
{
    public:
        NumaInfoTest() :
            BaseTest("numainfo_test")
        {
        }

        virtual void run() const
        {
#if defined linux
            unsigned nnodes = intern::num_nodes();

            TEST_CHECK(nnodes >= 1);

            unsigned lpus(sysconf(_SC_NPROCESSORS_CONF));
            unsigned * res = intern::cpu_to_node_array(nnodes, lpus);

            if (nnodes == 0)
            {
                TEST_CHECK(false);
            }
            else if (nnodes == 1)
            {
                TEST_CHECK(res == NULL);
            }
            else
            {
                TEST_CHECK(res != NULL);

                for (unsigned i(0) ; i < lpus ; ++i)
                {
                    TEST_CHECK(res[i] < lpus);
                }

                delete[] res;
            }

#else
            TEST_CHECK(true);
#endif
        }
} numainfo_test;
