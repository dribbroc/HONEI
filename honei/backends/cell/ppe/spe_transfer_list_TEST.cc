/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/backends/cell/ppe/spe_transfer_list.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace honei;
using namespace tests;

class SPETransferListTest :
    public QuickTest
{
    public:
        SPETransferListTest() :
            QuickTest("spe_transfer_list_test")
        {
        }

        virtual void run() const
        {
            float data[6] = { 0.0f, 0.5f, 2.3f, 4.6f, -1.2f, 3.5f };

            SPETransferList list(10, 10);

            for (unsigned i(0) ; i < 10 ; ++i)
            {
                list.add(data, 0);
                TEST_CHECK_EQUAL(list.size(), i + 1);
            }

            // Check if list size exceeding causes retval 0.
            SPETransferList::ListElement * retval = list.add(data, 0);
            TEST_CHECK_EQUAL(reinterpret_cast<unsigned long long>(retval), 0);

            SPETransferList list2(6, 20);

            TEST_CHECK_THROWS(list2.add(data, 6 * sizeof(float)), TransferListTransferSizeExceeded);

            // Check if two different EAHs cause retval 0.
            void * address(reinterpret_cast<void *>(0xFFFFFFFF00000000));
            void * address2(reinterpret_cast<void *>(0x0FFFFFFF00000000));
            SPETransferList list3(2, 1);
            list3.add(address, 0);
            SPETransferList::ListElement * retval3 = list3.add(address2, 0);
            TEST_CHECK_EQUAL(reinterpret_cast<unsigned long long>(retval3), 0);

        }
} spe_transfer_list_test;

class SPETransferListElementsTest :
    public QuickTest
{
    public:
        SPETransferListElementsTest() :
            QuickTest("spe_transfer_list_elements_test")
        {
        }

        virtual void run() const
        {
            char data[10];
            SPETransferList list(10, 10);

            SPETransferList::ListElement * j(list.elements());
            for (char * i(data), * i_end(data + 10) ; i != i_end ; ++i, ++j)
            {
                list.add(i, 1);
                TEST_CHECK_EQUAL(j->stall_and_notify, false);
                TEST_CHECK_EQUAL(j->reserved, 0);
                TEST_CHECK_EQUAL(j->transfer_size, 1);
                TEST_CHECK_EQUAL(j->effective_address_low, reinterpret_cast<unsigned long long>(i) & 0xFFFFFFFF);
            }
        }
} spe_transfer_list_elements_test;

