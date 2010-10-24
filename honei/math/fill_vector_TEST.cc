/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/fill_vector.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <endian_swap.hh>

//#include <cstdio>
//#include <cstdlib>

#include <fstream>
using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class FillVectorTest:
    public BaseTest
{
    private:
        unsigned long _size;
    public:
        FillVectorTest(const std::string & tag, unsigned long size) :
            BaseTest("FillVectorTestDirichlet1Neumann <" + tag + "> with size: " + stringify(size))
        {
            this->_size = size;
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> result(_size);
            FillVector<Tag_, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(result);
            result.lock(lm_read_only);
            std::cout << result << std::endl;
            TEST_CHECK(true);
            result.unlock(lm_read_only);
        }
};
FillVectorTest<tags::CPU, float> fill_matrix_float_2("float", 25ul);
FillVectorTest<tags::CPU, double> fill_matrix_double_2("double", 25ul);
