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
#include <unittest/unittest.hh>
#include <honei/lbm/source_grid.hh>
#include <honei/lbm/grid.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class SourceGridLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        SourceGridLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("source_grid_labswe_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            PackedGridData<lbm_lattice_types::D2Q9, DataType_> data;

            data.h = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            data.u = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            data.v = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            DataType_ g(9.81);
            DataType_ e(1.);
            DenseVector<DataType_> db(1000ul, DataType_(0.));
            DenseVector<DataType_> result(1000ul);

            SourceGrid<Tag_, lbm_applications::LABSWE, lbm_source_types::SIMPLE, lbm_source_schemes::BASIC>::
                value(data, result, db, g);
            for(unsigned long j(0); j < 1000; ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result[j], - g * 1.23456 * 0., std::numeric_limits<DataType_>::epsilon());
            }
        }
};
SourceGridLABSWETest<tags::CPU, float> source_grid_test_float("float");
SourceGridLABSWETest<tags::CPU, double> source_grid_test_double("double");
SourceGridLABSWETest<tags::CPU::MultiCore, float> source_grid_test_float_mc("float");
SourceGridLABSWETest<tags::CPU::MultiCore, double> source_grid_test_double_mc("double");
#ifdef HONEI_SSE
SourceGridLABSWETest<tags::CPU::SSE, float> source_grid_test_float_sse("float");
SourceGridLABSWETest<tags::CPU::SSE, double> source_grid_test_double_sse("double");
SourceGridLABSWETest<tags::CPU::MultiCore::SSE, float> source_grid_test_float_mc_sse("float");
SourceGridLABSWETest<tags::CPU::MultiCore::SSE, double> source_grid_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
SourceGridLABSWETest<tags::Cell, float> source_grid_test_float_cell("float");
SourceGridLABSWETest<tags::Cell, double> source_grid_test_double_cell("double");
#endif

