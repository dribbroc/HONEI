/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/dense_vector.hh>
#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/backends/mpi/operations.hh>
#include <honei/util/unittest.hh>

#include <string>
#include <limits>
#include <cmath>
#include <iostream>


using namespace honei;
using namespace tests;


template <typename DataType_>
class DenseVectorMPIQuickTest :
    public QuickTest
{
    public:
        DenseVectorMPIQuickTest(const std::string & type) :
            QuickTest("dense_vector_mpi_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseVector<DataType_> src(4711);
            for (unsigned long i(0) ; i < src.size() ; ++i)
                src[i] = i + 10;

            DenseVectorMPI<DataType_> dm0(src);
            DenseVectorMPI<DataType_> dm1(dm0);
            DenseVectorMPI<DataType_> dm2(dm1.copy());

            TEST_CHECK((unsigned long) dm0.vector().elements() == (unsigned long) dm1.vector().elements());
            TEST_CHECK((unsigned long) dm1.vector().elements() != (unsigned long) dm2.vector().elements());
            for (unsigned long i(0) ; i < dm2.local_size() ; ++i)
                TEST_CHECK_EQUAL(dm2[i], src[i + dm2.offset()]);
        }
};
DenseVectorMPIQuickTest<double> dense_vector_mpi_quick_test_double("double");

template <typename DataType_>
class DenseVectorMPIPartQuickTest :
    public QuickTest
{
    public:
        DenseVectorMPIPartQuickTest(const std::string & type) :
            QuickTest("dense_vector_mpi_part_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseVector<DataType_> src(4711);
            for (unsigned long i(0) ; i < src.size() ; ++i)
                src[i] = i + 10 * mpi::mpi_comm_rank();

            DenseVectorMPI<DataType_> dm0(src, mpi::mpi_comm_size() * src.size());
            DenseVectorMPI<DataType_> dm1(dm0);
            DenseVectorMPI<DataType_> dm2(dm1.copy());
            TEST_CHECK_EQUAL(dm2, dm0);

            for (unsigned long i(0) ; i < dm2.local_size() ; ++i)
                TEST_CHECK_EQUAL(dm2[i], src[i]);
        }
};
DenseVectorMPIPartQuickTest<double> dense_vector_mpi_part_quick_test_double("double");
