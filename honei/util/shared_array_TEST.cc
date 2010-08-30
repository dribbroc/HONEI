/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/util/shared_array.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using namespace tests;

template <typename DT_>
class SharedArrayAlignmentTest :
    public QuickTest
{
    public:
        SharedArrayAlignmentTest(const std::string & type) :
            QuickTest("shared_array_alignment_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            const unsigned repetitions(1000);
            const unsigned size(123);

            for (unsigned i(0) ; i < repetitions ; ++i)
            {
                SharedArray<float> array(size);
                unsigned long lsbs(reinterpret_cast<unsigned long>(array.get()));

                lsbs &= 0xf;
                TEST_CHECK_EQUAL(lsbs, 0ul);
            }
        }
};

SharedArrayAlignmentTest<float> shared_array_alignment_float_test("float");
SharedArrayAlignmentTest<double> shared_array_alignment_double_test("double");
