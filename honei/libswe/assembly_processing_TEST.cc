/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/libswe/assembly_processing.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <iostream>

using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class AssemblyProcessingTest:
    public BaseTest
{
    public:
        AssemblyProcessingTest(const std::string & tag) :
            BaseTest("AssemblyProcessing test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {

            DT1_ delta_t(2);
            DT1_ delta_x(2);
            DT1_ delta_y(2);
            unsigned long d_width = 1;
            unsigned long d_height = 1;

            unsigned long entries(3 * ((d_width * d_height) + 4 * (d_width + d_height + 4)));
            DenseVector<DT1_> u(entries, DT1_(5));
            DenseVector<DT1_> v(entries, DT1_(5));
            DenseVector<DT1_> w(entries, DT1_(5));
            DenseVector<DT1_> c(3, DT1_(10));
            DenseVector<DT1_> d(3, DT1_(10));

            BandedMatrix<DT1_> m1(entries);
            BandedMatrix<DT1_> m2(entries);
            BandedMatrix<DT1_> m3(entries);
            BandedMatrix<DT1_> m4(entries);
            BandedMatrix<DT1_> m6(entries);
            BandedMatrix<DT1_> m8(entries);


            AssemblyProcessing<tags::CPU, assembly_types::MAIN::M1M3>::value(m1, m3, u, v, delta_t, delta_x, d_width, d_height, c );
            AssemblyProcessing<tags::CPU, assembly_types::MAIN::M2M4>::value(m2, m4, u, w, delta_t, delta_y, d_width, d_height, d );
            AssemblyProcessing<tags::CPU, assembly_types::QUICK::M6>::value(m1, m6, c, d_width, d_height);
            AssemblyProcessing<tags::CPU, assembly_types::QUICK::M8>::value(m2, m8, d, d_width, d_height);

            cout << "Tested by visual verification M. Geveler 2007." << endl;
            TEST_CHECK(true);
        }
};
AssemblyProcessingTest<tags::CPU, float> assembly_test_float("float");
AssemblyProcessingTest<tags::CPU, double> assembly_test_double("double");

