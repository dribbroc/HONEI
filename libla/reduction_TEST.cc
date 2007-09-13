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
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DT_> * dv1 (new DenseVector<DT_>(size, DT_(2)));
                BandedMatrix<DT_> bm1(size, dv1);
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
            BandedMatrix<DT_> bm1(size, dv1);
            DenseVector<DT_> sum(Reduction<rt_sum>::value(bm1));

            TEST_CHECK_EQUAL(sum, *dv1);
        }
};
BandedMatrixReductionQuickTest<float> banded_matrix_reduction_to_sum_quick_test_float("float");
BandedMatrixReductionQuickTest<double> banded_matrix_reduction_to_sum_quick_test_double("double");

template <typename DT_>
class DenseMatrixReductionToSumTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToSumTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size + 1, DT_(1));
                DenseVector<DT_> dv1(size + 1, DT_(size + 1));
                DenseVector<DT_> sum(Reduction<rt_sum>::value(dm1));

                TEST_CHECK_EQUAL(sum, dv1);
            }
        }
};
DenseMatrixReductionToSumTest<float> dense_matrix_reduction_to_sum_test_float("float");
DenseMatrixReductionToSumTest<double> dense_matrix_reduction_to_sum_test_double("double");

template <typename DT_>
class DenseMatrixReductionQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DT_> dm1(size, size + 1, DT_(1));
            DenseVector<DT_> dv1(size + 1, DT_(size));
            DenseVector<DT_> sum(Reduction<rt_sum>::value(dm1));

            TEST_CHECK_EQUAL(sum, dv1);
        }
};
DenseMatrixReductionQuickTest<float> dense_matrix_reduction_to_sum_quick_test_float("float");
DenseMatrixReductionQuickTest<double> dense_matrix_reduction_to_sum_quick_test_double("double");

template <typename DT_>
class SparseMatrixReductionToSumTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToSumTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size + 1, size / 8 + 1);
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
                DenseVector<DT_> sum(Reduction<rt_sum>::value(sm1));

                TEST_CHECK_EQUAL(sum, dv1);
            }
        }
};
SparseMatrixReductionToSumTest<float> sparse_matrix_reduction_to_sum_test_float("float");
SparseMatrixReductionToSumTest<double> sparse_matrix_reduction_to_sum_test_double("double");

template <typename DT_>
class SparseMatrixReductionQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DT_> sm1(size, size + 1, size / 8 + 1);
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
            DenseVector<DT_> sum(Reduction<rt_sum>::value(sm1));

            TEST_CHECK_EQUAL(sum, dv1);
        }
};
SparseMatrixReductionQuickTest<float> sparse_matrix_reduction_to_sum_quick_test_float("float");
SparseMatrixReductionQuickTest<double> sparse_matrix_reduction_to_sum_quick_test_double("double");

template <typename DT_>
class DenseVectorReductionToSumTest :
    public BaseTest
{
    public:
        DenseVectorReductionToSumTest(const std::string & type) :
            BaseTest("dense_vector_eduction_to_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DT_> dv(size);
                for (typename Vector<DT_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
                }

                DT_ v1(Reduction<rt_sum>::value(dv));
                DT_ s1(size * (size + 1) / 2 / 1.23456789);
                // Behavious similar to size^2 * eps
                DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

DenseVectorReductionToSumTest<float> dense_vector_reduction_to_sum_test_float("float");
DenseVectorReductionToSumTest<double> dense_vector_reduction_to_sum_test_double("double");

template <typename DT_>
class DenseVectorReductionToSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DT_> dv(size);
            for (typename Vector<DT_>::ElementIterator i(dv.begin_elements()), i_end(dv.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DT_>((i.index() + 1) / 1.23456789);
            }

            DT_ v1(Reduction<rt_sum>::value(dv));
            DT_ s1(size * (size + 1) / 2 / 1.23456789);
            // Behavious similar to size^2 * eps
            DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
        }
};
DenseVectorReductionToSumQuickTest<float>  dense_vector_reduction_to_sum_quick_test_float("float");
DenseVectorReductionToSumQuickTest<double> dense_vector_reduction_to_sum_quick_test_double("double");

template <typename DT_>
class SparseVectorReductionToSumTest :
    public BaseTest
{
    public:
        SparseVectorReductionToSumTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DT_ s1(0);
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = static_cast<DT_>((i.index() +1) / 1.23456789);
                    s1 += *i;
                }

                DT_ v1(Reduction<rt_sum>::value(sv1));
                DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

SparseVectorReductionToSumTest<float> sparse_vector_reduction_to_sum_test_float("float");
SparseVectorReductionToSumTest<double> sparse_vector_reduction_to_sum_test_double("double");

template <typename DT_>
class SparseVectorReductionToSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            DT_ s1(0);
            SparseVector<DT_> sv1(size, size / 8 + 1);
            for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = static_cast<DT_>((i.index() +1) / 1.23456789);
                s1 += *i;
            }

            DT_ v1(Reduction<rt_sum>::value(sv1));
            DT_ eps1(s1 * 10 * std::numeric_limits<DT_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
        }
};

SparseVectorReductionToSumQuickTest<float> sparse_reduction_to_sum_quick_test_float("float");
SparseVectorReductionToSumQuickTest<double> sparse_reduction_to_sum_quick_test_double("double");

template <typename DT_>
class DenseMatrixReductionToMinTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToMinTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_min_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_min>::value(dm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
            }
        }
};

DenseMatrixReductionToMinTest<float> dense_matrix_reduction_to_min_test_float("float");
DenseMatrixReductionToMinTest<double> dense_matrix_reduction_to_min_test_double("double");

template <typename DT_>
class DenseMatrixReductionToMinQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionToMinQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_min_quick_test<" + type + ">")
        {
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

            DenseVector<DT_> v1(Reduction<rt_min>::value(dm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
        }
};

DenseMatrixReductionToMinQuickTest<float> dense_matrix_reduction_to_min_quick_test_float("float");
DenseMatrixReductionToMinQuickTest<double> dense_matrix_reduction_to_min_quick_test_double("double");

template <typename DT_>
class SparseMatrixReductionToMinTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToMinTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_min_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_min>::value(sm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
            }
        }
};

SparseMatrixReductionToMinTest<float> sparse_matrix_reduction_to_min_test_float("float");
SparseMatrixReductionToMinTest<double> sparse_matrix_reduction_to_min_test_double("double");

template <typename DT_>
class SparseMatrixReductionToMinQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionToMinQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_min_quick_test<" + type + ">")
        {
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

            DenseVector<DT_> v1(Reduction<rt_min>::value(sm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size - size));
        }
};

SparseMatrixReductionToMinQuickTest<float> sparse_matrix_reduction_to_min_quick_test_float("float");
SparseMatrixReductionToMinQuickTest<double> sparse_matrix_reduction_to_min_quick_test_double("double");

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
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
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

template <typename DT_>
class DenseVectorReductionToMinTest :
    public BaseTest
{
    public:
        DenseVectorReductionToMinTest(const std::string & type) :
            BaseTest("dense_vector_reduction_to_min_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size);
                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_min>::value(dv1));
                TEST_CHECK_EQUAL(v1, 0);
            }
        }
};

DenseVectorReductionToMinTest<float> dense_vector_reduction_to_min_test_float("float");
DenseVectorReductionToMinTest<double> dense_vector_reduction_to_min_test_double("double");

template <typename DT_>
class DenseVectorReductionToMinQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToMinQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_min_quick_test<" + type + ">")
        {
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

            DT_ v1(Reduction<rt_min>::value(dv1));
            TEST_CHECK_EQUAL(v1, 0);
        }
};

DenseVectorReductionToMinQuickTest<float> dense_reduction_to_min_quick_test_float("float");
DenseVectorReductionToMinQuickTest<double> dense_reduction_to_min_quick_test_double("double");

template <typename DT_>
class SparseVectorReductionToMinTest :
    public BaseTest
{
    public:
        SparseVectorReductionToMinTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_min_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_min>::value(sv1));
                TEST_CHECK_EQUAL(v1, 0);
            }
        }
};

SparseVectorReductionToMinTest<float> sparse_vector_reduction_to_min_test_float("float");
SparseVectorReductionToMinTest<double> sparse_vector_reduction_to_min_test_double("double");

template <typename DT_>
class SparseVectorReductionToMinQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToMinQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_min_quick_test<" + type + ">")
        {
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

            DT_ v1(Reduction<rt_min>::value(sv1));
            TEST_CHECK_EQUAL(v1, 0);
        }
};

SparseVectorReductionToMinQuickTest<float> sparse_reduction_to_min_quick_test_float("float");
SparseVectorReductionToMinQuickTest<double> sparse_reduction_to_min_quick_test_double("double");

template <typename DT_>
class DenseMatrixReductionToMaxTest :
    public BaseTest
{
    public:
        DenseMatrixReductionToMaxTest(const std::string & type) :
            BaseTest("dense_matrix_reduction_to_max_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseMatrix<DT_> dm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(dm1.begin_elements()), i_end(dm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_max>::value(dm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
            }
        }
};

DenseMatrixReductionToMaxTest<float> dense_matrix_reduction_to_max_test_float("float");
DenseMatrixReductionToMaxTest<double> dense_matrix_reduction_to_max_test_double("double");

template <typename DT_>
class DenseMatrixReductionToMaxQuickTest :
    public QuickTest
{
    public:
        DenseMatrixReductionToMaxQuickTest(const std::string & type) :
            QuickTest("dense_matrix_reduction_to_max_quick_test<" + type + ">")
        {
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

            DenseVector<DT_> v1(Reduction<rt_max>::value(dm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
        }
};

DenseMatrixReductionToMaxQuickTest<float> dense_matrix_reduction_to_max_quick_test_float("float");
DenseMatrixReductionToMaxQuickTest<double> dense_matrix_reduction_to_max_quick_test_double("double");

template <typename DT_>
class SparseMatrixReductionToMaxTest :
    public BaseTest
{
    public:
        SparseMatrixReductionToMaxTest(const std::string & type) :
            BaseTest("sparse_matrix_reduction_to_max_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseMatrix<DT_> sm1(size, size);
                for (typename MutableMatrix<DT_>::ElementIterator i(sm1.begin_elements()), i_end(sm1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DenseVector<DT_> v1(Reduction<rt_max>::value(sm1));
                TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
            }
        }
};

SparseMatrixReductionToMaxTest<float> sparse_matrix_reduction_to_max_test_float("float");
SparseMatrixReductionToMaxTest<double> sparse_matrix_reduction_to_max_test_double("double");

template <typename DT_>
class SparseMatrixReductionToMaxQuickTest :
    public QuickTest
{
    public:
        SparseMatrixReductionToMaxQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_reduction_to_max_quick_test<" + type + ">")
        {
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

            DenseVector<DT_> v1(Reduction<rt_max>::value(sm1));
            TEST_CHECK_EQUAL(v1[size-1], (size*size)-1);
        }
};

SparseMatrixReductionToMaxQuickTest<float> sparse_matrix_reduction_to_max_quick_test_float("float");
SparseMatrixReductionToMaxQuickTest<double> sparse_matrix_reduction_to_max_quick_test_double("double");

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
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
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

template <typename DT_>
class DenseVectorReductionToMaxTest :
    public BaseTest
{
    public:
        DenseVectorReductionToMaxTest(const std::string & type) :
            BaseTest("dense_vector_reduction_to_max_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DT_> dv1(size);
                for (typename Vector<DT_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_max>::value(dv1));
                TEST_CHECK_EQUAL(v1, size-1);
            }
        }
};

DenseVectorReductionToMaxTest<float> dense_vector_reduction_to_max_test_float("float");
DenseVectorReductionToMaxTest<double> dense_vector_reduction_to_max_test_double("double");

template <typename DT_>
class DenseVectorReductionToMaxQuickTest :
    public QuickTest
{
    public:
        DenseVectorReductionToMaxQuickTest(const std::string & type) :
            QuickTest("dense_vector_reduction_to_max_quick_test<" + type + ">")
        {
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

            DT_ v1(Reduction<rt_max>::value(dv1));
            TEST_CHECK_EQUAL(v1, size-1);
        }
};

DenseVectorReductionToMaxQuickTest<float> dense_reduction_to_max_quick_test_float("float");
DenseVectorReductionToMaxQuickTest<double> dense_reduction_to_max_quick_test_double("double");


template <typename DT_>
class SparseVectorReductionToMaxTest :
    public BaseTest
{
    public:
        SparseVectorReductionToMaxTest(const std::string & type) :
            BaseTest("sparse_vector_reduction_to_max_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DT_> sv1(size, size / 8 + 1);
                for (typename Vector<DT_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = i.index();
                }

                DT_ v1(Reduction<rt_max>::value(sv1));
                TEST_CHECK_EQUAL(v1, size-1);
            }
        }
};

SparseVectorReductionToMaxTest<float> sparse_vector_reduction_to_max_test_float("float");
SparseVectorReductionToMaxTest<double> sparse_vector_reduction_to_max_test_double("double");

template <typename DT_>
class SparseVectorReductionToMaxQuickTest :
    public QuickTest
{
    public:
        SparseVectorReductionToMaxQuickTest(const std::string & type) :
            QuickTest("sparse_vector_reduction_to_max_quick_test<" + type + ">")
        {
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

            DT_ v1(Reduction<rt_max>::value(sv1));
            TEST_CHECK_EQUAL(v1, size-1);
        }
};

SparseVectorReductionToMaxQuickTest<float> sparse_reduction_to_max_quick_test_float("float");
SparseVectorReductionToMaxQuickTest<double> sparse_reduction_to_max_quick_test_double("double");

