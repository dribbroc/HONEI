/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@tu-dortmund.de>
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

#include <honei/libmath/hessenberg.hh>
#include <honei/libmath/qr_decomposition.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DT_>
class QRDecompositionTest :
    public QuickTest
{
    private:

    public:
        QRDecompositionTest(const std::string & type) :
            QuickTest("qr_decomposition_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(6);

            DenseMatrix<DT_> dm(size, size, DT_(1));

            Hessenberg<Tag_>::value(dm);
            QRDecomposition<Tag_>::value(dm);

            // No checks yet.
        }
};

QRDecompositionTest<tags::CPU, float> qr_decomposition_test_float("float");

