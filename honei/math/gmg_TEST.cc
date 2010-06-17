/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.tu-dortmund.de>
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

#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <cmath>
#include <honei/math/gmg.hh>
#include <honei/math/problem_factory.hh>
#include <honei/math/methods.hh>
#include <honei/math/preconditioning.hh>
#include <honei/math/multigrid.hh>
#include <honei/la/sum.hh>
#include <honei/util/operation_wrapper.hh>

using namespace honei;
using namespace tests;
using namespace std;


template <typename Tag_, typename DT1_>
class GMGTest:
    public BaseTest
{
    public:

        static void factory(std::vector<SparseMatrixELL<DT1_> > & a, std::vector<SparseMatrixELL<DT1_> > & b, std::vector<SparseMatrixELL<DT1_> > & c, std::vector<DenseVector<DT1_> > & d, unsigned long min, unsigned long max)
        {
            DenseVector<DT1_> blub(10, 5);
            d.push_back(blub);
        }

        GMGTest(const std::string & tag) :
        BaseTest("GMG test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::vector<unsigned long> test_cycle;
            FileFactory<DT1_, SparseMatrixELL<DT1_>, SparseMatrixELL<DT1_> > factory("bla");
            GMGInfo<DT1_, SparseMatrixELL<DT1_>, SparseMatrixELL<DT1_> > info(GMGInfoFactory<DT1_, SparseMatrixELL<DT1_>, SparseMatrixELL<DT1_>, Jacobi<Tag_>, Preconditioning<Tag_, methods::NONE>, ConjugateGradients<Tag_, methods::NONE>, Prolongation<Tag_, methods::NONE>, Restriction<Tag_, methods::NONE>, Norm<vnt_l_two, false, Tag_> >::create(factory.factory, 1ul, 1ul, 1ul, 1ul, 1ul, 1ul, 1ul, DT1_(0), DT1_(0), DT1_(0), test_cycle));
            GMG<Tag_, methods::NONE>::value(info);
            std::cout<<info.rhs_vectors.front();

            /*std::vector<std::tr1::function<void ()> > functors;
            DenseVector<DT1_> bla(100, DT1_(4711));
            DenseVector<DT1_> blup(100, DT1_(4713));

            OperationWrapper<Sum<Tag_>, DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT1_> > sum_operation(blup);
            std::tr1::function<void ()> func(std::tr1::bind(sum_operation, bla, blup));
            functors.push_back(func);

            functors.front()();
            std::cout << bla;*/

        }
};
GMGTest<tags::CPU, float> cpu_gmg_test_float("float");

