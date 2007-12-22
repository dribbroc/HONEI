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
                        j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 20 == 0)
                    {
                        *i = DataType_(-0.25) * (i.index() + 1);
                        *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                    }
                    else if (i.index() % 30 == 0)
                    {
                        *i = DataType_(0);
                        *j = DataType_(0);
                    }
                    else if (i.index() % 33 == 0)
                    {
                        *i = DataType_(-0);
                        *j = DataType_(-0);
                    }
                    else
                    {
                        *i = (i.index() + 1) / DataType_(1.234);
                        *j = 1 / ((i.index() + 1) / DataType_(1.234));
                    }
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
#ifdef HONEI_SSE
DenseVectorElementInverseTest<tags::CPU::SSE, float> sse_dense_vector_element_inverse_test_float("SSE float");
DenseVectorElementInverseTest<tags::CPU::SSE, double> sse_dense_vector_element_inverse_test_double("SSE double");
DenseVectorElementInverseTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_element_inverse_test_float("MC SSE float");
DenseVectorElementInverseTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_element_inverse_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseVectorElementInverseTest<tags::Cell, float> cell_dense_vector_element_inverse_test_float("Cell float");
#endif

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
            unsigned long size(4711);
            DenseVector<DataType_> dv1(size, DataType_(0)), dv2(size, DataType_(0));
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 20 == 0)
                {
                    *i = DataType_(-0.25) * (i.index() + 1);
                    *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                }
                else if (i.index() % 30 == 0)
                {
                    *i = DataType_(0);
                    *j = DataType_(0);
                }
                else if (i.index() % 33 == 0)
                {
                    *i = DataType_(-0);
                    *j = DataType_(-0);
                }
                else
                {
                    *i = (i.index() + 1) / DataType_(1.234);
                    *j = 1 / ((i.index() + 1) / DataType_(1.234));
                }
            }
            ElementInverse<Tag_>::value(dv1);
            TEST_CHECK_EQUAL(dv1, dv2);
        }
};
DenseVectorElementInverseQuickTest<tags::CPU, float> dense_vector_element_inverse_quick_test_float("float");
DenseVectorElementInverseQuickTest<tags::CPU, double> dense_vector_element_inverse_quick_test_double("double");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_element_inverse_quick_test_float("MC float");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_element_inverse_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorElementInverseQuickTest<tags::CPU::SSE, float> sse_dense_vector_element_inverse_quick_test_float("SSE float");
DenseVectorElementInverseQuickTest<tags::CPU::SSE, double> sse_dense_vector_element_inverse_quick_test_double("SSE double");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_element_inverse_quick_test_float("MC SSE float");
DenseVectorElementInverseQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_element_inverse_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseVectorElementInverseQuickTest<tags::Cell, float> cell_dense_vector_element_inverse_quick_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeElementInverseTest :
    public BaseTest
{
    public:
        DenseVectorRangeElementInverseTest(const std::string & type) :
            BaseTest("dense_vector_range_element_inverse_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                for (int o(0) ; o < 4 ; ++o)
                {
                    DenseVector<DataType_> dv1(size * 2, DataType_(0)), dv2(size, DataType_(0));
                    DenseVectorRange<DataType_> dv1r(dv1, size, o + ((size -1) / 4));
                    for (typename Vector<DataType_>::ElementIterator i(dv1r.begin_elements()), i_end(dv1r.end_elements()),
                            j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                    {
                        if (i.index() % 20 == 0)
                        {
                            *i = DataType_(-0.25) * (i.index() + 1);
                            *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                        }
                        else if (i.index() % 30 == 0)
                        {
                            *i = DataType_(0);
                            *j = DataType_(0);
                        }
                        else if (i.index() % 33 == 0)
                        {
                            *i = DataType_(-0);
                            *j = DataType_(-0);
                        }
                        else
                        {
                            *i = (i.index() + 1) / DataType_(1.234);
                            *j = 1 / ((i.index() + 1) / DataType_(1.234));
                        }
                    }
                    ElementInverse<Tag_>::value(dv1r);
                    TEST_CHECK_EQUAL(dv1r, dv2);
                }
            }
        }
};

DenseVectorRangeElementInverseTest<tags::CPU, float> dense_vector_range_element_inverse_test_float("float");
DenseVectorRangeElementInverseTest<tags::CPU, double> dense_vector_range_element_inverse_test_double("double");
DenseVectorRangeElementInverseTest<tags::CPU::MultiCore, float> mc_dense_vector_range_element_inverse_test_float("MC float");
DenseVectorRangeElementInverseTest<tags::CPU::MultiCore, double> mc_dense_vector_range_element_inverse_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeElementInverseTest<tags::CPU::SSE, float> sse_dense_vector_range_element_inverse_test_float("SSE float");
DenseVectorRangeElementInverseTest<tags::CPU::SSE, double> sse_dense_vector_range_element_inverse_test_double("SSE double");
DenseVectorRangeElementInverseTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_element_inverse_test_float("MC SSE float");
DenseVectorRangeElementInverseTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_element_inverse_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
//DenseVectorRangeElementInverseTest<tags::Cell, float> cell_dense_vector_range_element_inverse_test_float("Cell float");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeElementInverseQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeElementInverseQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_element_inverse_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(4711);
            for (int o(0) ; o < 4 ; ++o)
            {
                DenseVector<DataType_> dv1(size * 2, DataType_(0)), dv2(size, DataType_(0));
                DenseVectorRange<DataType_> dv1r(dv1, size, o + ((size -1) / 4));
                for (typename Vector<DataType_>::ElementIterator i(dv1r.begin_elements()), i_end(dv1r.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 20 == 0)
                    {
                        *i = DataType_(-0.25) * (i.index() + 1);
                        *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                    }
                    else if (i.index() % 30 == 0)
                    {
                        *i = DataType_(0);
                        *j = DataType_(0);
                    }
                    else if (i.index() % 33 == 0)
                    {
                        *i = DataType_(-0);
                        *j = DataType_(-0);
                    }
                    else
                    {
                        *i = (i.index() + 1) / DataType_(1.234);
                        *j = 1 / ((i.index() + 1) / DataType_(1.234));
                    }
                }
                ElementInverse<Tag_>::value(dv1r);
                TEST_CHECK_EQUAL(dv1r, dv2);
            }
        }
};

DenseVectorRangeElementInverseQuickTest<tags::CPU, float> dense_vector_range_element_inverse_quick_test_float("float");
DenseVectorRangeElementInverseQuickTest<tags::CPU, double> dense_vector_range_element_inverse_quick_test_double("double");
DenseVectorRangeElementInverseQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_range_element_inverse_quick_test_float("MC float");
DenseVectorRangeElementInverseQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_element_inverse_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeElementInverseQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_element_inverse_quick_test_float("SSE float");
DenseVectorRangeElementInverseQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_element_inverse_quick_test_double("SSE double");
DenseVectorRangeElementInverseQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_vector_range_element_inverse_quick_test_float("MC SSE float");
DenseVectorRangeElementInverseQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_vector_range_element_inverse_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
//DenseVectorRangeElementInverseQuickTest<tags::Cell, float> cell_dense_vector_range_element_inverse_quick_test_float("Cell float");
#endif

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
                sv1[0] = DataType_(0);
                sv2[0] = DataType_(0);

                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sv1), sv2);
            }
        }
};
SparseVectorElementInverseTest<tags::CPU, float> sparse_vector_element_inverse_test_float("float");
SparseVectorElementInverseTest<tags::CPU, double> sparse_vector_element_inverse_test_double("double");
SparseVectorElementInverseTest<tags::CPU::MultiCore, float> mc_sparse_vector_element_inverse_test_float("MC float");
SparseVectorElementInverseTest<tags::CPU::MultiCore, double> mc_sparse_vector_element_inverse_test_double("MC double");
#ifdef HONEI_SSE
SparseVectorElementInverseTest<tags::CPU::SSE, float> sse_sparse_vector_element_inverse_test_float("SSE float");
SparseVectorElementInverseTest<tags::CPU::SSE, double> sse_sparse_vector_element_inverse_test_double("SSE double");
SparseVectorElementInverseTest<tags::CPU::MultiCore::SSE, float> sse_mc_sparse_vector_element_inverse_test_float("MC SSE float");
SparseVectorElementInverseTest<tags::CPU::MultiCore::SSE, double> sse_mc_sparse_vector_element_inverse_test_double("MC SSE double");
#endif

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
            sv1[0] = DataType_(0);
            sv2[0] = DataType_(0);

            TEST_CHECK_EQUAL(ElementInverse<>::value(sv1), sv2);
        }
};
SparseVectorElementInverseQuickTest<tags::CPU, float> sparse_vector_element_inverse_quick_test_float("float");
SparseVectorElementInverseQuickTest<tags::CPU, double> sparse_vector_element_inverse_quick_test_double("double");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore, float> mc_sparse_vector_element_inverse_quick_test_float("MC float");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore, double> mc_sparse_vector_element_inverse_quick_test_double("MC double");
#ifdef HONEI_SSE
SparseVectorElementInverseQuickTest<tags::CPU::SSE, float> sse_sparse_vector_element_inverse_quick_test_float("SSE float");
SparseVectorElementInverseQuickTest<tags::CPU::SSE, double> sse_sparse_vector_element_inverse_quick_test_double("SSE double");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_sparse_vector_element_inverse_quick_test_float("MC SSE float");
SparseVectorElementInverseQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_sparse_vector_element_inverse_quick_test_double("MC SSE double");
#endif

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
                    *i = i.index() - 1;
                    if (i.index() - 1 == 0) *j = DataType_(0);
                    else *j = 1 / DataType_(i.index() - 1);
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
                    *i = i.index() - 1;
                    if (i.index() - 1 == 0) *j = DataType_(0);
                    else *j = 1 / DataType_(i.index() - 1);
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
                DenseMatrix<DataType_> dm1(size, size + 3, DataType_(0)),
                        dm2(size, size + 3, DataType_(0));
                for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()),
                        j(dm2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 20 == 0)
                    {
                        *i = DataType_(-0.25) * (i.index() + 1);
                        *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                    }
                    else if (i.index() % 30 == 0)
                    {
                        *i = DataType_(0);
                        *j = DataType_(0);
                    }
                    else if (i.index() % 33 == 0)
                    {
                        *i = DataType_(-0);
                        *j = DataType_(-0);
                    }
                    else
                    {
                        *i = (i.index() + 1) / DataType_(1.234);
                        *j = 1 / ((i.index() + 1) / DataType_(1.234));
                    }
                }

                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(dm1), dm2);
            }
        }
};
DenseMatrixElementInverseTest<tags::CPU, float> dense_matrix_element_inverse_test_float("float");
DenseMatrixElementInverseTest<tags::CPU, double> dense_matrix_element_inverse_test_double("double");
DenseMatrixElementInverseTest<tags::CPU::MultiCore, float> mc_dense_matrix_element_inverse_test_float("MC float");
DenseMatrixElementInverseTest<tags::CPU::MultiCore, double> mc_dense_matrix_element_inverse_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixElementInverseTest<tags::CPU::SSE, float> sse_dense_matrix_element_inverse_test_float("SSE float");
DenseMatrixElementInverseTest<tags::CPU::SSE, double>sse_dense_matrix_element_inverse_test_double("SSE double");
DenseMatrixElementInverseTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_element_inverse_test_float("MC SSE float");
DenseMatrixElementInverseTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_element_inverse_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixElementInverseTest<tags::Cell, float> cell_dense_matrix_element_inverse_test_float("Cell float");
#endif

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
            unsigned long size(333);
            DenseMatrix<DataType_> dm1(size, size + 3, DataType_(0)),
                    dm2(size, size + 3, DataType_(0));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()),
                    j(dm2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 20 == 0)
                {
                    *i = DataType_(-0.25) * (i.index() + 1);
                    *j = 1 / ((DataType_(-0.25) * (i.index() + 1)));
                }
                else if (i.index() % 30 == 0)
                {
                    *i = DataType_(0);
                    *j = DataType_(0);
                }
                else if (i.index() % 33 == 0)
                {
                    *i = DataType_(-0);
                    *j = DataType_(-0);
                }
                else
                {
                    *i = (i.index() + 1) / DataType_(1.234);
                    *j = 1 / ((i.index() + 1) / DataType_(1.234));
                }
            }

            TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(dm1), dm2);
        }
};
DenseMatrixElementInverseQuickTest<tags::CPU, float>  dense_matrix_element_inverse_quick_test_float("float");
DenseMatrixElementInverseQuickTest<tags::CPU, double> dense_matrix_element_inverse_quick_test_double("double");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_element_inverse_quick_test_float("MC float");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_element_inverse_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseMatrixElementInverseQuickTest<tags::CPU::SSE, float>  sse_dense_matrix_element_inverse_quick_test_float("SSE float");
DenseMatrixElementInverseQuickTest<tags::CPU::SSE, double> sse_dense_matrix_element_inverse_quick_test_double("SSE double");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_dense_matrix_element_inverse_quick_test_float("MC SSE float");
DenseMatrixElementInverseQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_dense_matrix_element_inverse_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CELL
DenseMatrixElementInverseQuickTest<tags::Cell, float> cell_dense_matrix_element_inverse_quick_test_float("Cell float");
#endif

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
                sm1[0][0] = DataType_(0);
                sm2[0][0] = DataType_(0);
                TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sm1), sm2);
            }
        }
};
SparseMatrixElementInverseTest<tags::CPU, float> sparse_matrix_element_inverse_test_float("float");
SparseMatrixElementInverseTest<tags::CPU, double> sparse_matrix_element_inverse_test_double("double");
SparseMatrixElementInverseTest<tags::CPU::MultiCore, float> mc_sparse_matrix_element_inverse_test_float("MC float");
SparseMatrixElementInverseTest<tags::CPU::MultiCore, double> mc_sparse_matrix_element_inverse_test_double("MC double");
#ifdef HONEI
SparseMatrixElementInverseTest<tags::CPU::SSE, float> sse_sparse_matrix_element_inverse_test_float("SSE float");
SparseMatrixElementInverseTest<tags::CPU::SSE, double> sse_sparse_matrix_element_inverse_test_double("SSE double");
SparseMatrixElementInverseTest<tags::CPU::MultiCore::SSE, float> sse_mc_sparse_matrix_element_inverse_test_float("MC SSE float");
SparseMatrixElementInverseTest<tags::CPU::MultiCore::SSE, double> sse_mc_sparse_matrix_element_inverse_test_double("MC SSE double");
#endif

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
            sm1[0][0] = DataType_(0);
            sm2[0][0] = DataType_(0);
            TEST_CHECK_EQUAL(ElementInverse<Tag_>::value(sm1), sm2);
        }
};
SparseMatrixElementInverseQuickTest<tags::CPU, float> sparse_matrix_element_inverse_quick_test_float("float");
SparseMatrixElementInverseQuickTest<tags::CPU, double> sparse_matrix_element_inverse_quick_test_double("double");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_element_inverse_quick_test_float("MC float");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_element_inverse_quick_test_double("MC double");
#ifdef HONEI_SSE
SparseMatrixElementInverseQuickTest<tags::CPU::SSE, float> sse_sparse_matrix_element_inverse_quick_test_float("SSE float");
SparseMatrixElementInverseQuickTest<tags::CPU::SSE, double> sse_sparse_matrix_element_inverse_quick_test_double("SSE double");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore::SSE, float> sse_mc_sparse_matrix_element_inverse_quick_test_float("MC SSE float");
SparseMatrixElementInverseQuickTest<tags::CPU::MultiCore::SSE, double> sse_mc_sparse_matrix_element_inverse_quick_test_double("MC SSE double");
#endif
