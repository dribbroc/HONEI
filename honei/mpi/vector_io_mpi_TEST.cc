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
#include <honei/math/vector_io.hh>
#include <honei/mpi/vector_io_mpi.hh>

#include <string>
#include <limits>
#include <cmath>
#include <iostream>


using namespace honei;
using namespace tests;


template <typename DataType_>
class DenseVectorMPIIOQuickTest :
    public QuickTest
{
    public:
        DenseVectorMPIIOQuickTest(const std::string & type) :
            QuickTest("dense_vector_mpi_io_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/poisson_advanced4/sort_0/rhs_4";
            DenseVector<DataType_> srhs(VectorIO<io_formats::EXP>::read_vector(filename, DataType_(0)));
            DenseVectorMPI<DataType_> ref_rhs(srhs);
            DenseVectorMPI<DataType_> rhs(VectorIOMPI<io_formats::EXP>::read_vector(filename, DataType_(0)));

            TEST_CHECK_EQUAL(rhs, ref_rhs);
        }
};
DenseVectorMPIIOQuickTest<double> dense_vector_mpi_io_quick_test_double("double");
