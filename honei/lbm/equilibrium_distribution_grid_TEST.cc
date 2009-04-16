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
#include <honei/lbm/equilibrium_distribution_grid.hh>

#include <limits>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class EqDisGridLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        EqDisGridLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("eq_dis_grid_labswe_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            PackedGridData<lbm_lattice_types::D2Q9, DataType_> data;
            PackedGridInfo<lbm_lattice_types::D2Q9> info;

            data.h = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            info.limits = new DenseVector<unsigned long>(2);
            unsigned long begin(0);
            unsigned long end(data.h->size());
            (*info.limits)[0] = begin;
            (*info.limits)[1] = end;
            data.u = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            data.v = new DenseVector<DataType_>(1000ul, DataType_(1.23456));
            DataType_ g(9.81);
            DataType_ e(1.);
            data.distribution_x = new DenseVector<DataType_>(9ul, DataType_(2.));
            data.distribution_y = new DenseVector<DataType_>(9ul, DataType_(2.));
            data.f_eq_0 = new DenseVector<DataType_>(1000);
            data.f_eq_1 = new DenseVector<DataType_>(1000);
            data.f_eq_2 = new DenseVector<DataType_>(1000);
            data.f_eq_3 = new DenseVector<DataType_>(1000);
            data.f_eq_4 = new DenseVector<DataType_>(1000);
            data.f_eq_5 = new DenseVector<DataType_>(1000);
            data.f_eq_6 = new DenseVector<DataType_>(1000);
            data.f_eq_7 = new DenseVector<DataType_>(1000);
            data.f_eq_8 = new DenseVector<DataType_>(1000);
            EquilibriumDistributionGrid<Tag_, lbm_applications::LABSWE>::
                value(g, e, info, data);

            data.f_eq_0->lock(lm_read_only);
            data.f_eq_1->lock(lm_read_only);
            data.f_eq_2->lock(lm_read_only);
            data.f_eq_3->lock(lm_read_only);
            data.f_eq_4->lock(lm_read_only);
            data.f_eq_5->lock(lm_read_only);
            data.f_eq_6->lock(lm_read_only);
            data.f_eq_7->lock(lm_read_only);
            data.f_eq_8->lock(lm_read_only);
            for(unsigned long i(begin); i < end; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_eq_0)[i], (1.23456 - ((5. * 9.81 * 1.23456 * 1.23456) / (6.)) - ((2. * 1.23456) / (3.)) * ((1.23456 * 1.23456) + (1.23456 * 1.23456))), std::numeric_limits<DataType_>::epsilon() * 50.);

                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_eq_1)[i], ((9.81 * 1.23456 * 1.23456) / 6. + ((1.23456 / 3.) * (2. * 1.23456 + 2. * 1.23456)) + ((1.23456 / 2.) * (2. * 1.23456 * 2. * 1.23456 + 2. * 2. * 1.23456 * 2. * 1.23456 + 2. * 1.23456 * 2. * 1.23456)) - ((1.23456 / 6.) * (1.23456 * 1.23456 + 1.23456 * 1.23456))), std::numeric_limits<DataType_>::epsilon() * 50.);
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_eq_2)[i], ((9.81 * 1.23456 * 1.23456) / 24. + ((1.23456 / 12.) * (2. * 1.23456 + 2. * 1.23456)) + ((1.23456 / 8.) * (2. * 1.23456 * 2. * 1.23456 + 2. * 2. * 1.23456 * 2. * 1.23456 + 2. * 1.23456 * 2. * 1.23456)) - ((1.23456 / 24.) * (1.23456 * 1.23456 + 1.23456 * 1.23456))), std::numeric_limits<DataType_>::epsilon() * 50.);
            }

            data.f_eq_0->unlock(lm_read_only);
            data.f_eq_1->unlock(lm_read_only);
            data.f_eq_2->unlock(lm_read_only);
            data.f_eq_3->unlock(lm_read_only);
            data.f_eq_4->unlock(lm_read_only);
            data.f_eq_5->unlock(lm_read_only);
            data.f_eq_6->unlock(lm_read_only);
            data.f_eq_7->unlock(lm_read_only);
            data.f_eq_8->unlock(lm_read_only);
        }
};
EqDisGridLABSWETest<tags::CPU, float> eq_dist_grid_test_float("float");
EqDisGridLABSWETest<tags::CPU, double> eq_dist_grid_test_double("double");
#ifdef HONEI_SSE
EqDisGridLABSWETest<tags::CPU::SSE, float> sse_eq_dist_grid_test_float("float");
EqDisGridLABSWETest<tags::CPU::SSE, double> sse_eq_dist_grid_test_double("double");
#endif
#ifdef HONEI_CUDA
EqDisGridLABSWETest<tags::GPU::CUDA, float> cuda_eq_dist_grid_test_float("float");
#endif
#ifdef HONEI_CELL
EqDisGridLABSWETest<tags::Cell, float> cell_eq_dist_grid_test_float("float");
#endif
