/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. Math is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * Math is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/math/transposition.hh>
#include <unittest/unittest.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class TranspositionTEST:
    public BaseTest
{
    public:
        TranspositionTEST(const std::string & tag) :
            BaseTest("Transposition Test<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            SparseMatrix<DT1_> source(5, 10);
            SparseMatrix<DT1_> target(10, 5);
            source(1,9) = DT1_(1);
            source(4,2) = DT1_(2);
            source(2,8) = DT1_(3);

            SparseMatrix<DT1_>ref(10 , 5);
            ref(9, 1) = DT1_(1);
            ref(2, 4) = DT1_(2);
            ref(8, 2) = DT1_(3);

            Transposition<Tag_>::value(source, target);
            TEST_CHECK_EQUAL(target, ref);
        }
};

