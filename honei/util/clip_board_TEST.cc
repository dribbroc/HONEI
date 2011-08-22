/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/unittest.hh>
#include <honei/util/clip_board.hh>

using namespace honei;
using namespace tests;

class ClipBoardQuickTest :
    public QuickTest
{
    public:
        ClipBoardQuickTest() :
            QuickTest("clip_board_quick_test")
        {
        }

        virtual void run() const
        {
            for (double i(0) ; i < 5 ; ++i)
                ClipBoard<double, 1>::instance()->push_back(i);

            for (double i(4) ; i >=0 ; --i)
                TEST_CHECK_EQUAL((ClipBoard<double, 1>::instance()->at(i)), i);

            for (double i(0) ; i < 5 ; ++i)
                ClipBoard<double, 1>::instance()->pop_back();

            TEST_CHECK_EQUAL((ClipBoard<double, 1>::instance()->size()), 0ul);
        }
} clip_board_quick_test;
