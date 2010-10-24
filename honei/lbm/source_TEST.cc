/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/lbm/tags.hh>
#include <honei/util/unittest.hh>
#include <honei/lbm/source.hh>

#include <limits>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class SourceLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        SourceLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("source_labswe_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> h(1000ul, 1000ul, DataType_(1.23456));
            DenseMatrix<DataType_> db(1000ul, 1000ul, DataType_(0.));
            DataType_ g(9.81);
            DenseMatrix<DataType_> result(1000ul, 1000ul);

            Source<Tag_, lbm_applications::LABSWE, lbm_force::SIMPLE, lbm_source_schemes::BASIC>::
                value(result, h, db, g);
            for(unsigned long i(0); i < 1000; ++i)
            {
                for(unsigned long j(0); j < 1000; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(result(i,j), - g * 1.23456 * 0., std::numeric_limits<DataType_>::epsilon());
                }
            }
            Source<Tag_, lbm_applications::LABSWE, lbm_force::CONSTANT, lbm_source_schemes::BASIC>::
                value(result, DataType_(0.000024));
            for(unsigned long i(0); i < 1000; ++i)
            {
                for(unsigned long j(0); j < 1000; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(result(i,j), 0.000024, std::numeric_limits<DataType_>::epsilon());
                }
            }
        }
};
SourceLABSWETest<tags::CPU, float> source_test_float("float");
SourceLABSWETest<tags::CPU, double> source_test_double("double");
SourceLABSWETest<tags::CPU::MultiCore, float> source_test_float_mc("float");
SourceLABSWETest<tags::CPU::MultiCore, double> source_test_double_mc("double");
#ifdef HONEI_SSE
SourceLABSWETest<tags::CPU::SSE, float> source_test_float_sse("float");
SourceLABSWETest<tags::CPU::SSE, double> source_test_double_sse("double");
SourceLABSWETest<tags::CPU::MultiCore::SSE, float> source_test_float_mc_sse("float");
SourceLABSWETest<tags::CPU::MultiCore::SSE, double> source_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
SourceLABSWETest<tags::Cell, float> source_test_float_cell("float");
SourceLABSWETest<tags::Cell, double> source_test_double_cell("double");
#endif

