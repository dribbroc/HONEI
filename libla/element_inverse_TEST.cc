/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <libla/element_inverse.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace honei;
using  namespace tests;

template <typename Tag_, typename DataType_>
class DenseVectorElementInverseTest :
    public BaseTest
{
    public:
        DenseVectorElementInverseTest(const std::string & type) :
            BaseTest("dense_vector_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0));
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                    ++j;
                }
                ElementInverse<Tag_>::value(dv1);
                TEST_CHECK_EQUAL(dv1, dv2);
            }
        }
};

DenseVectorElementInverseTest<tags::CPU, float> dense_vector_element_inverse_test_float("float");
DenseVectorElementInverseTest<tags::CPU, double> dense_vector_element_inverse_test_double("double");
DenseVectorElementInverseTest<tags::CPU::MultiCore, float> mc_dense_vector_element_inverse_test_float("MC float");
DenseVectorElementInverseTest<tags::CPU::MultiCore, double> mc_dense_vector_element_inverse_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseVectorElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_vector_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0));
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()) ; i != i_end ; ++i)
            {
                *i = i.index() + 1;
                *j = 1 / static_cast<DataType_>(i.index() + 1);
                ++j;
            }
            ElementInverse<Tag_>::value(dv1);
            TEST_CHECK_EQUAL(dv1, dv2);
        }
};
DenseVectorElementInverseQuickTest<tags::CPU, float> dense_vector_element_inverse_quick_test_float("float");
DenseVectorElementInverseQuickTest<tags::CPU, double> dense_vector_element_inverse_quick_test_double("double");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_element_inverse_quick_test_float("MC float");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_element_inverse_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseVectorElementInverseTest :
    public BaseTest
{
    public:
        SparseVectorElementInverseTest(const std::string & type) :
            BaseTest("sparse_vector_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                        j(sv2.begin_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = i.index() + 1;
                        *j = 1 / static_cast<DataType_>(i.index() + 1);
                    }
                    ++j;
                }

                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sv1), sv2);
            }
        }
};
SparseVectorElementInverseTest<tags::CPU, float> sparse_vector_element_inverse_test_float("float");
SparseVectorElementInverseTest<tags::CPU, double> sparse_vector_element_inverse_test_double("double");
SparseVectorElementInverseTest<tags::CPU::MultiCore, float> mc_sparse_vector_element_inverse_test_float("MC float");
SparseVectorElementInverseTest<tags::CPU::MultiCore, double> mc_sparse_vector_element_inverse_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseVectorElementInverseQuickTest :
    public QuickTest
{
    public:
        SparseVectorElementInverseQuickTest(const std::string & type) :
            QuickTest("sparse_vector_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (20);
            SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = i.index() + 1;
                    *j = 1 / static_cast<DataType_>(i.index() + 1);
                }
                ++j;
            }

            TEST_CHECK_EQUAL(ElementInverse<>::value(sv1), sv2);
        }
};
SparseVectorElementInverseQuickTest<tags::CPU, float> sparse_vector_element_inverse_quick_test_float("float");
SparseVectorElementInverseQuickTest<tags::CPU, double> sparse_vector_element_inverse_quick_test_double("double");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore, float> mc_sparse_vector_element_inverse_quick_test_float("MC float");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore, double> mc_sparse_vector_element_inverse_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class BandedMatrixElementInverseTest :
    public BaseTest
{
    public:
        BandedMatrixElementInverseTest(const std::string & type) :
            BaseTest("banded_matrix_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size);
                DenseVector<DataType_> dv2(size);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / DataType_(i.index() + 1);
                    ++j;
                }
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(bm1), bm2);
            }
        }
};
BandedMatrixElementInverseTest<tags::CPU, float> banded_matrix_element_inverse_test_float("float");
BandedMatrixElementInverseTest<tags::CPU, double> banded_matrix_element_inverse_test_double("double");
BandedMatrixElementInverseTest<tags::CPU::MultiCore, float> mc_banded_matrix_element_inverse_test_float("MC float");
BandedMatrixElementInverseTest<tags::CPU::MultiCore, double> mc_banded_matrix_element_inverse_test_double("MC double");

template <typename Tag_, typename DataType_>
class BandedMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("banded_matrix_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
                unsigned long size(11);
                DenseVector<DataType_> dv1(size);
                DenseVector<DataType_> dv2(size);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / DataType_(i.index() + 1);
                    ++j;
                }
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);

                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(bm1), bm2);
        }
};
BandedMatrixElementInverseQuickTest<tags::CPU, float> banded_matrix_element_inverse_quick_test_float("float");
BandedMatrixElementInverseQuickTest<tags::CPU, double> banded_matrix_element_inverse_quick_test_double("double");
BandedMatrixElementInverseQuickTest<tags::CPU::MultiCore, float> mc_banded_matrix_element_inverse_quick_test_float("MC float");
BandedMatrixElementInverseQuickTest<tags::CPU::MultiCore, double> mc_banded_matrix_element_inverse_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseMatrixElementInverseTest :
    public BaseTest
{
    public:
        DenseMatrixElementInverseTest(const std::string & type) :
            BaseTest("dense_matrix_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 8) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size, DataType_(0)),
                        dm2(size, size, DataType_(0));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()),
                        j(dm2.begin_elements()) ; i != i_end ; ++i)
                {
                    *i = i.index() + 1;
                    *j = 1 / DataType_(i.index() + 1);
                    ++j;
                }

                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(dm1), dm2);
            }
        }
};
DenseMatrixElementInverseTest<tags::CPU, float> dense_matrix_element_inverse_test_float("float");
DenseMatrixElementInverseTest<tags::CPU, double> dense_matrix_element_inverse_test_double("double");
DenseMatrixElementInverseTest<tags::CPU::MultiCore, float> mc_dense_matrix_element_inverse_test_float("MC float");
DenseMatrixElementInverseTest<tags::CPU::MultiCore, double> mc_dense_matrix_element_inverse_test_double("MC double");

template <typename Tag_, typename DataType_>
class DenseMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_matrix_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> dm1(3, 2, DataType_(0)),
                    dm2(3, 2, DataType_(0));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()),
                    j(dm2.begin_elements()) ; i != i_end ; ++i)
            {
                *i = i.index() + 1;
                *j = 1 / static_cast<DataType_>(i.index() + 1);
                ++j;
            }

            TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(dm1), dm2);
        }
};
DenseMatrixElementInverseQuickTest<tags::CPU, float>  dense_matrix_element_inverse_quick_test_float("float");
DenseMatrixElementInverseQuickTest<tags::CPU, double> dense_matrix_element_inverse_quick_test_double("double");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_element_inverse_quick_test_float("MC float");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_element_inverse_quick_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseMatrixElementInverseTest :
    public BaseTest
{
    public:
        SparseMatrixElementInverseTest(const std::string & type) :
            BaseTest("sparse_matrix_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 8) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                    sm2(size+1, size, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), j(sm2.begin_elements());
                    i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = i.index() + 1;
                        *j = 1 / DataType_(i.index() + 1);

                    }
                }
                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sm1), sm2);
            }
        }
};
SparseMatrixElementInverseTest<tags::CPU, float> sparse_matrix_element_inverse_test_float("float");
SparseMatrixElementInverseTest<tags::CPU, double> sparse_matrix_element_inverse_test_double("double");
SparseMatrixElementInverseTest<tags::CPU::MultiCore, float> mc_sparse_matrix_element_inverse_test_float("MC float");
SparseMatrixElementInverseTest<tags::CPU::MultiCore, double> mc_sparse_matrix_element_inverse_test_double("MC double");

template <typename Tag_, typename DataType_>
class SparseMatrixElementInverseQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementInverseQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (7);
            SparseMatrix<DataType_> sm1(size+1, size, size / 8 + 1),
                sm2(size+1, size, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), j(sm2.begin_elements());
                i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0)
                {
                    *i = i.index() + 1;
                    *j = 1 / DataType_(i.index() + 1);

                }
            }
            TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sm1), sm2);
        }
};
SparseMatrixElementInverseQuickTest<tags::CPU, float> sparse_matrix_element_inverse_quick_test_float("float");
SparseMatrixElementInverseQuickTest<tags::CPU, double> sparse_matrix_element_inverse_quick_test_double("double");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_element_inverse_quick_test_float("MC float");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_element_inverse_quick_test_double("MC double");
