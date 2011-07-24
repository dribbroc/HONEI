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
            DenseVector<DataType_> src(4711, DataType_(7));
            DenseVectorMPI<DataType_> dm00(src, 0, 5);
            DenseVectorMPI<DataType_> dm0(dm00);
            DenseVectorMPI<DataType_> dm11(src, 1, 5);
            DenseVectorMPI<DataType_> dm1(dm11.copy());
            DenseVectorMPI<DataType_> dm2(src, 2, 5);
            DenseVectorMPI<DataType_> dm3(src, 3, 5);
            DenseVectorMPI<DataType_> dm4(src, 4, 5);
            std::cout<<dm0.size()<<" "<<dm0.offset()<<std::endl;
            std::cout<<dm1.size()<<" "<<dm1.offset()<<std::endl;
            std::cout<<dm2.size()<<" "<<dm2.offset()<<std::endl;
            std::cout<<dm3.size()<<" "<<dm3.offset()<<std::endl;
            std::cout<<dm4.size()<<" "<<dm4.offset()<<std::endl;
        }
};
DenseVectorMPIQuickTest<float>  dense_vector_mpi_quick_test_float("float");
DenseVectorMPIQuickTest<double> dense_vector_mpi_quick_test_double("double");
