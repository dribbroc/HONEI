/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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


#include <honei/math/spai.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>


using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT_>
class SpaiTestSparse:
    public BaseTest
{
    public:
        SpaiTestSparse(const std::string & tag) :
            BaseTest("Spai test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            SparseMatrix<DT_> sm(10, 10);
            for (unsigned long i(0) ; i < 10 ; i++)
            {
                sm(i,i) = 3;
            }
            sm(1, 0) = 5.123;
            sm(9, 8) = -5.123;
            sm(8, 0) = 1.123;
            std::cout<<sm;
            SparseMatrix<DT_> m(SPAI::value(sm));
            std::cout<<m;

        }
};

SpaiTestSparse<tags::CPU, float> spai_test_sparse_ell_float("float");
SpaiTestSparse<tags::CPU, double> spai_test_sparse_ell_double("double");
