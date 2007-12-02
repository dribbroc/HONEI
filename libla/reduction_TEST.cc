/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libla/dense_matrix.hh>
#include <libla/reduction.hh>
#include <matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using namespace tests;

template <typename DT_>
class BandedMatrixReductionToSumTest :
    public BaseTest
{
    public:
        BandedMatrixReductionToSumTest(const std::string & type) :
            BaseTest("banded_matrix_reduction_to_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseVector<DT_> * dv1 (new DenseVector<DT_>(size, DT_(2)));
                BandedMatrix<DT_> bm1(size, *dv1);
                DenseVector<DT_> sum(Reduction<rt_sum>::value(bm1));

                TEST_CHECK_EQUAL(sum, *dv1);
            }
        }
};
BandedMatrixReductionToSumTest<float> banded_matrix_reduction_to_sum_test_float("float");
BandedMatrixReductionToSumTest<double> banded_matrix_reduction_to_sum_test_double("double");

template <typename DT_>
class BandedMatrixReductionQuickTest :
    public QuickTest
{
    public:
        BandedMatrixReductionQuickTest(const std::string & type) :
            QuickTest("banded_matrix_reduction_to_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DT_> * dv1 (new DenseVector<DT_>(size, DT_(2)));
            BandedMatrix<DT_> bm1(size, *dv1);
            DenseVector<DT_> sum(Reduction<rt_sum>::value(bm1));

            TEST_CHECK_EQUAL(sum, *dv1);
        }
};
BandedMatrixReductionQuickTest<float> banded_matrix_reduction_to_sum_quick_test_float("float");
BandedMatrixReductionQuickTest<double> banded_matrix_reduction_to_sum_quick_test_double("double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionToSumTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToSumTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 9) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size+1, size, DT_(1));
                DenseVector<DT_> dv1(size + 1, DT_(size ));
                DenseVector<DT_> sum(Reduction<rt_sum, Tag_>::value(dm1));

                TEST_CHECK_EQUAL(sum, dv1);
            }
        }
};
DenseMatrixReductionToSumTest<tags::CPU, float> dense_matrix_reduction_to_sum_test_float("float");
DenseMatrixReductionToSumTest<tags::CPU, double> dense_matrix_reduction_to_sum_test_double("double");
DenseMatrixReductionToSumTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_sum_test_float("MC float");
DenseMatrixReductionToSumTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_sum_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DT_> dm1(size+1, size, DT_(1));
            DenseVector<DT_> dv1(size + 1, DT_(size));
            DenseVector<DT_> sum(Reduction<rt_sum, Tag_>::value(dm1));

            TEST_CHECK_EQUAL(sum, dv1);
        }
};
DenseMatrixReductionQuickTest<tags::CPU, float> dense_matrix_reduction_to_sum_quick_test_float("float");
DenseMatrixReductionQuickTest<tags::CPU, double> dense_matrix_reduction_to_sum_quick_test_double("double");
DenseMatrixReductionQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_sum_quick_test_float("MC float");
DenseMatrixReductionQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_sum_quick_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseMatrixReductionToSumTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToSumTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 9) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size+1, size, size / 8 + 1);
                DenseVector<DT_> dv1(size + 1, DT_(0));
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = DT_(2);
                        dv1[i.row()] += DT_(2);
                    }
                }
                DenseVector<DT_> sum(Reduction<rt_sum, Tag_>::value(sm1));

                TEST_CHECK_EQUAL(sum, dv1);
            }
        }
};
SparseMatrixReductionToSumTest<tags::CPU, float> sparse_matrix_reduction_to_sum_test_float("float");
SparseMatrixReductionToSumTest<tags::CPU, double> sparse_matrix_reduction_to_sum_test_double("double");
SparseMatrixReductionToSumTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_sum_test_float("MC float");
SparseMatrixReductionToSumTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_sum_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseMatrixReductionQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DT_> sm1(size+1, size, size / 8 + 1);
            DenseVector<DT_> dv1(size + 1, DT_(0));
            for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = DT_(2);
                    dv1[i.row()] += DT_(2);
                }
            }
            DenseVector<DT_> sum(Reduction<rt_sum, Tag_>::value(sm1));

            TEST_CHECK_EQUAL(sum, dv1);
        }
};
SparseMatrixReductionQuickTest<tags::CPU, float> sparse_matrix_reduction_to_sum_quick_test_float("float");
SparseMatrixReductionQuickTest<tags::CPU, double> sparse_matrix_reduction_to_sum_quick_test_double("double");
SparseMatrixReductionQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_sum_quick_test_float("MC float");
SparseMatrixReductionQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_sum_quick_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseVectorReductionToSumTest :
    public BaseTest
{
    public:
        DenseVectorReductionToSumTest(const std::string & type) :
            BaseTest("dense_vector_reduction_to_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DT_> dv(size);
                for (typename Vector<DT_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
                }

                DT_ v1(Reduction<rt_sum, Tag_>::value(dv));
                DT_ s1(size * (size + 1) / 2 / 1.23456789);
                // Behaviour similar to size^2 * eps
                DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

DenseVectorReductionToSumTest<tags::CPU, float> dense_vector_reduction_to_sum_test_float("float");
DenseVectorReductionToSumTest<tags::CPU, double> dense_vector_reduction_to_sum_test_double("double");
DenseVectorReductionToSumTest<tags::CPU::MultiCore, float> mc_dense_vector_reduction_to_sum_test_float("MC float");
DenseVectorReductionToSumTest<tags::CPU::MultiCore, double> mc_dense_vector_reduction_to_sum_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseVectorReductionToSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DT_> dv(size);
            for (typename Vector<DT_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
            }

            DT_ v1(Reduction<rt_sum, Tag_>::value(dv));
            DT_ s1(size * (size + 1) / 2 / 1.23456789);
            // Behaviour similar to size^2 * eps
            DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
        }
};
DenseVectorReductionToSumQuickTest<tags::CPU, float>  dense_vector_reduction_to_sum_quick_test_float("float");
DenseVectorReductionToSumQuickTest<tags::CPU, double> dense_vector_reduction_to_sum_quick_test_double("double");
DenseVectorReductionToSumQuickTest<tags::CPU::MultiCore, float>  mc_dense_vector_reduction_to_sum_quick_test_float("MC float");
DenseVectorReductionToSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_reduction_to_sum_quick_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorReductionToSumQuickTest<tags::Cell, float> dense_vector_reduction_to_sum_quick_test_float_cell("Cell float");
#endif


template <typename Tag_, typename DT_>
class SparseVectorReductionToSumTest :
    public BaseTest
{
    public:
        SparseVectorReductionToSumTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DT_ s1(0);
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0)
                    {
                        *i = static_cast<DT_>((i.index() +1) / 1.23456789);
                        s1 += *i;
                    }                
                }

                DT_ v1(Reduction<rt_sum, Tag_>::value(sv1));
                DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

SparseVectorReductionToSumTest<tags::CPU, float> sparse_vector_reduction_to_sum_test_float("float");
SparseVectorReductionToSumTest<tags::CPU, double> sparse_vector_reduction_to_sum_test_double("double");
SparseVectorReductionToSumTest<tags::CPU::MultiCore, float> mc_sparse_vector_reduction_to_sum_test_float("MC float");
SparseVectorReductionToSumTest<tags::CPU::MultiCore, double> mc_sparse_vector_reduction_to_sum_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseVectorReductionToSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DT_ s1(0);
            SparseVector<DT_> sv1(size, size / 8 + 1);
            for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0)
                {
                    *i = static_cast<DT_>((i.index() +1) / 1.23456789);
                    s1 += *i;
                }
            }

            DT_ v1(Reduction<rt_sum, Tag_>::value(sv1));
            DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
        }
};

SparseVectorReductionToSumQuickTest<tags::CPU, float> sparse_reduction_to_sum_quick_test_float("float");
SparseVectorReductionToSumQuickTest<tags::CPU, double> sparse_reduction_to_sum_quick_test_double("double");
SparseVectorReductionToSumQuickTest<tags::CPU::MultiCore, float> mc_sparse_reduction_to_sum_quick_test_float("MC float");
SparseVectorReductionToSumQuickTest<tags::CPU::MultiCore, double> mc_sparse_reduction_to_sum_quick_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionToMinTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToMinTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_min_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_min,Tag_>::value(dm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
            }
        }
};
DenseMatrixReductionToMinTest<tags::CPU, float> dense_matrix_reduction_to_min_test_float("float");
DenseMatrixReductionToMinTest<tags::CPU, double> dense_matrix_reduction_to_min_test_double("double");
DenseMatrixReductionToMinTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_min_test_float("MC float");
DenseMatrixReductionToMinTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_min_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionToMinQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionToMinQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_min_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseMatrix<DT_> dm1(size, size);
            for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_min, Tag_>::value(dm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
        }
};
DenseMatrixReductionToMinQuickTest<tags::CPU, float> dense_matrix_reduction_to_min_quick_test_float("float");
DenseMatrixReductionToMinQuickTest<tags::CPU, double> dense_matrix_reduction_to_min_quick_test_double("double");
DenseMatrixReductionToMinQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_min_quick_test_float("MC float");
DenseMatrixReductionToMinQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_min_quick_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseMatrixReductionToMinTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToMinTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_min_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_min, Tag_>::value(sm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
            }
        }
};
SparseMatrixReductionToMinTest<tags::CPU, float> sparse_matrix_reduction_to_min_test_float("float");
SparseMatrixReductionToMinTest<tags::CPU, double> sparse_matrix_reduction_to_min_test_double("double");
SparseMatrixReductionToMinTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_min_test_float("MC float");
SparseMatrixReductionToMinTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_min_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseMatrixReductionToMinQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionToMinQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_min_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DT_> sm1(size, size);
            for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_min, Tag_>::value(sm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
        }
};
SparseMatrixReductionToMinQuickTest<tags::CPU, float> sparse_matrix_reduction_to_min_quick_test_float("float");
SparseMatrixReductionToMinQuickTest<tags::CPU, double> sparse_matrix_reduction_to_min_quick_test_double("double");
SparseMatrixReductionToMinQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_min_quick_test_float("MC float");
SparseMatrixReductionToMinQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_min_quick_test_double("MC double");

template <typename DT_>
class BandedMatrixReductionToMinTest :
    public BaseTest
{
    public:
        BandedMatrixReductionToMinTest(const std::string & type) :
            BaseTest("banded_matrix_reduction_to_min_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                BandedMatrix<DT_> bm1(size);
                for (typename MutableMatrix<DT_>::ElementIterator i(bm1.begin_elements()), i_end(bm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }


                DenseVector<DT_> v1(Reduction<rt_min>::value(bm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
            }
        }
};

BandedMatrixReductionToMinTest<float> banded_matrix_reduction_to_min_test_float("float");
BandedMatrixReductionToMinTest<double> banded_matrix_reduction_to_min_test_double("double");

template <typename DT_>
class BandedMatrixReductionToMinQuickTest :
    public QuickTest
{
    public:
        BandedMatrixReductionToMinQuickTest(const std::string & type) :
            QuickTest("banded_matrix_reduction_to_min_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            BandedMatrix<DT_> bm1(size);
            for (typename MutableMatrix<DT_>::ElementIterator i(bm1.begin_elements()), i_end(bm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_min>::value(bm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
        }
};

BandedMatrixReductionToMinQuickTest<float> banded_matrix_reduction_to_min_quick_test_float("float");
BandedMatrixReductionToMinQuickTest<double> banded_matrix_reduction_to_min_quick_test_double("double");

template <typename Tag_, typename DT_>
class DenseVectorReductionToMinTest :
    public BaseTest
{
    public:
        DenseVectorReductionToMinTest(const std::string & type) :
            BaseTest("dense_vector_reduction_to_min_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size);
                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_min, Tag_>::value(dv1));
                TEST_CHECK_EQUAL(v1, 0);
            }
        }
};

DenseVectorReductionToMinTest<tags::CPU, float> dense_vector_reduction_to_min_test_float("float");
DenseVectorReductionToMinTest<tags::CPU, double> dense_vector_reduction_to_min_test_double("double");
DenseVectorReductionToMinTest<tags::CPU::MultiCore, float> mc_dense_vector_reduction_to_min_test_float("MC float");
DenseVectorReductionToMinTest<tags::CPU::MultiCore, double> mc_dense_vector_reduction_to_min_test_double("MC double");


template <typename Tag_, typename DT_>
class DenseVectorReductionToMinQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToMinQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_min_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DT_> dv1(size);
            for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DT_ v1(Reduction<rt_min, Tag_>::value(dv1));
            TEST_CHECK_EQUAL(v1, 0);
        }
};

DenseVectorReductionToMinQuickTest<tags::CPU, float> dense_reduction_to_min_quick_test_float("float");
DenseVectorReductionToMinQuickTest<tags::CPU, double> dense_reduction_to_min_quick_test_double("double");
DenseVectorReductionToMinQuickTest<tags::CPU::MultiCore, float> mc_dense_reduction_to_min_quick_test_float("MC float");
DenseVectorReductionToMinQuickTest<tags::CPU::MultiCore, double> mc_dense_reduction_to_min_quick_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorReductionToMinQuickTest<tags::Cell, float> dense_vector_reduction_to_min_quick_test_float_cell("Cell float");
#endif





template <typename Tag_, typename DT_>
class SparseVectorReductionToMinTest :
    public BaseTest
{
    public:
        SparseVectorReductionToMinTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_min_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_min, Tag_>::value(sv1));
                TEST_CHECK_EQUAL(v1, 0);
            }
        }
};

SparseVectorReductionToMinTest<tags::CPU, float> sparse_vector_reduction_to_min_test_float("float");
SparseVectorReductionToMinTest<tags::CPU, double> sparse_vector_reduction_to_min_test_double("double");
SparseVectorReductionToMinTest<tags::CPU::MultiCore, float> mc_sparse_vector_reduction_to_min_test_float("MC float");
SparseVectorReductionToMinTest<tags::CPU::MultiCore, double> mc_sparse_vector_reduction_to_min_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseVectorReductionToMinQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToMinQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_min_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseVector<DT_> sv1(size, size / 8 + 1);
            for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DT_ v1(Reduction<rt_min, Tag_>::value(sv1));
            TEST_CHECK_EQUAL(v1, 0);
        }
};

SparseVectorReductionToMinQuickTest<tags::CPU, float> sparse_reduction_to_min_quick_test_float("float");
SparseVectorReductionToMinQuickTest<tags::CPU, double> sparse_reduction_to_min_quick_test_double("double");
SparseVectorReductionToMinQuickTest<tags::CPU::MultiCore, float> mc_sparse_reduction_to_min_quick_test_float("MC float");
SparseVectorReductionToMinQuickTest<tags::CPU::MultiCore, double> mc_sparse_reduction_to_min_quick_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionToMaxTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToMaxTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_max_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_max, Tag_>::value(dm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
            }
        }
};

DenseMatrixReductionToMaxTest<tags::CPU, float> dense_matrix_reduction_to_max_test_float("float");
DenseMatrixReductionToMaxTest<tags::CPU, double> dense_matrix_reduction_to_max_test_double("double");
DenseMatrixReductionToMaxTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_max_test_float("MC float");
DenseMatrixReductionToMaxTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_max_test_double("MC double");

template <typename Tag_, typename DT_>
class DenseMatrixReductionToMaxQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionToMaxQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_max_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseMatrix<DT_> dm1(size, size);
            for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_max, Tag_>::value(dm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
        }
};

DenseMatrixReductionToMaxQuickTest<tags::CPU, float> dense_matrix_reduction_to_max_quick_test_float("float");
DenseMatrixReductionToMaxQuickTest<tags::CPU, double> dense_matrix_reduction_to_max_quick_test_double("double");
DenseMatrixReductionToMaxQuickTest<tags::CPU::MultiCore, float> mc_dense_matrix_reduction_to_max_quick_test_float("MC float");
DenseMatrixReductionToMaxQuickTest<tags::CPU::MultiCore, double> mc_dense_matrix_reduction_to_max_quick_test_double("MC double");


template <typename Tag_, typename DT_>
class SparseMatrixReductionToMaxTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToMaxTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_max_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_max, Tag_>::value(sm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
            }
        }
};

SparseMatrixReductionToMaxTest<tags::CPU, float> sparse_matrix_reduction_to_max_test_float("float");
SparseMatrixReductionToMaxTest<tags::CPU, double> sparse_matrix_reduction_to_max_test_double("double");
SparseMatrixReductionToMaxTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_max_test_float("MC float");
SparseMatrixReductionToMaxTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_max_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseMatrixReductionToMaxQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionToMaxQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_max_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DT_> sm1(size, size);
            for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_max, Tag_>::value(sm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
        }
};

SparseMatrixReductionToMaxQuickTest<tags::CPU, float> sparse_matrix_reduction_to_max_quick_test_float("float");
SparseMatrixReductionToMaxQuickTest<tags::CPU, double> sparse_matrix_reduction_to_max_quick_test_double("double");
SparseMatrixReductionToMaxQuickTest<tags::CPU::MultiCore, float> mc_sparse_matrix_reduction_to_max_quick_test_float("MC float");
SparseMatrixReductionToMaxQuickTest<tags::CPU::MultiCore, double> mc_sparse_matrix_reduction_to_max_quick_test_double("MC double");

template <typename DT_>
class BandedMatrixReductionToMaxTest :
    public BaseTest
{
    public:
        BandedMatrixReductionToMaxTest(const std::string & type) :
            BaseTest("banded_matrix_reduction_to_max_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                BandedMatrix<DT_> bm1(size);
                for (typename MutableMatrix<DT_>::ElementIterator i(bm1.begin_elements()), i_end(bm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_max>::value(bm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
            }
        }
};

BandedMatrixReductionToMaxTest<float> banded_matrix_reduction_to_max_test_float("float");
BandedMatrixReductionToMaxTest<double> banded_matrix_reduction_to_max_test_double("double");

template <typename DT_>
class BandedMatrixReductionToMaxQuickTest :
    public QuickTest
{
    public:
        BandedMatrixReductionToMaxQuickTest(const std::string & type) :
            QuickTest("banded_matrix_reduction_to_max_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            BandedMatrix<DT_> bm1(size);
            for (typename MutableMatrix<DT_>::ElementIterator i(bm1.begin_elements()), i_end(bm1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DenseVector<DT_> v1(Reduction<rt_max>::value(bm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
        }
};

BandedMatrixReductionToMaxQuickTest<float> banded_matrix_reduction_to_max_quick_test_float("float");
BandedMatrixReductionToMaxQuickTest<double> banded_matrix_reduction_to_max_quick_test_double("double");

template <typename Tag_, typename DT_>
class DenseVectorReductionToMaxTest :
    public BaseTest
{
    public:
        DenseVectorReductionToMaxTest(const std::string & type) :
            BaseTest("dense_vector_reduction_to_max_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size);
                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_max, Tag_>::value(dv1));
                TEST_CHECK_EQUAL(v1, size-1);
            }
        }
};

DenseVectorReductionToMaxTest<tags::CPU, float> dense_vector_reduction_to_max_test_float("float");
DenseVectorReductionToMaxTest<tags::CPU, double> dense_vector_reduction_to_max_test_double("double");
DenseVectorReductionToMaxTest<tags::CPU::MultiCore, float> mc_dense_vector_reduction_to_max_test_float("MC float");
DenseVectorReductionToMaxTest<tags::CPU::MultiCore, double> mc_dense_vector_reduction_to_max_test_double("MC double");


template <typename Tag_, typename DT_>
class DenseVectorReductionToMaxQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToMaxQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_max_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DT_> dv1(size);
            for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DT_ v1(Reduction<rt_max,Tag_>::value(dv1));
            TEST_CHECK_EQUAL(v1, size-1);
        }
};

DenseVectorReductionToMaxQuickTest<tags::CPU, float> dense_reduction_to_max_quick_test_float("float");
DenseVectorReductionToMaxQuickTest<tags::CPU, double> dense_reduction_to_max_quick_test_double("double");
DenseVectorReductionToMaxQuickTest<tags::CPU::MultiCore, float> mc_dense_reduction_to_max_quick_test_float("MC float");
DenseVectorReductionToMaxQuickTest<tags::CPU::MultiCore, double> mc_dense_reduction_to_max_quick_test_double("MC double");
#ifdef HONEI_CELL
DenseVectorReductionToMaxQuickTest<tags::Cell, float> dense_vector_reduction_to_max_quick_test_float_cell("Cell float");
#endif






template <typename Tag_, typename DT_>
class SparseVectorReductionToMaxTest :
    public BaseTest
{
    public:
        SparseVectorReductionToMaxTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_max_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_max, Tag_>::value(sv1));
                TEST_CHECK_EQUAL(v1, size-1);
            }
        }
};

SparseVectorReductionToMaxTest<tags::CPU, float> sparse_vector_reduction_to_max_test_float("float");
SparseVectorReductionToMaxTest<tags::CPU, double> sparse_vector_reduction_to_max_test_double("double");
SparseVectorReductionToMaxTest<tags::CPU::MultiCore, float> mc_sparse_vector_reduction_to_max_test_float("MC float");
SparseVectorReductionToMaxTest<tags::CPU::MultiCore, double> mc_sparse_vector_reduction_to_max_test_double("MC double");

template <typename Tag_, typename DT_>
class SparseVectorReductionToMaxQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToMaxQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_max_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseVector<DT_> sv1(size, size / 8 + 1);
            for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = i.index();
            }

            DT_ v1(Reduction<rt_max, Tag_>::value(sv1));
            TEST_CHECK_EQUAL(v1, size-1);
        }
};

SparseVectorReductionToMaxQuickTest<tags::CPU, float> sparse_reduction_to_max_quick_test_float("float");
SparseVectorReductionToMaxQuickTest<tags::CPU, double> sparse_reduction_to_max_quick_test_double("double");
SparseVectorReductionToMaxQuickTest<tags::CPU::MultiCore, float> mc_sparse_reduction_to_max_quick_test_float("MC float");
SparseVectorReductionToMaxQuickTest<tags::CPU::MultiCore, double> mc_sparse_reduction_to_max_quick_test_double("MC double");

