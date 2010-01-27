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

#include <honei/backends/multicore/topology.hh>
#include <unittest/unittest.hh>

using namespace honei::mc;
using namespace tests;

class TopologyTest :
    public BaseTest
{
    public:
        TopologyTest() :
            BaseTest("topology_test")
        {
        }

        virtual void run() const
        {
            Topology * t(new Topology);

#if defined linux
            TEST_CHECK(t->num_lpus() == (unsigned) sysconf(_SC_NPROCESSORS_CONF));
#else
            TEST_CHECK(t->num_lpus() > 0u);
#endif

#if defined(__i386__) || defined(__x86_64__)
            TEST_CHECK(t->num_cpus() > 0u);
            TEST_CHECK(t->num_cores() > 0u);
            TEST_CHECK(t->ht_factor() > 0u);
#endif
            TEST_CHECK(true);
        }
} topology_test;
