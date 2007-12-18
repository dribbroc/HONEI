/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBLA_GUARD_SUM_MC_HH
#define LIBLA_GUARD_SUM_MC_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libutil/tags.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>

namespace honei
{
    // Forward declaration.
    template <typename Tag_> struct Sum;

    template <typename Tag_> struct MCSum
    {
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & a, const SparseVector<DT2_> & b, unsigned long offset)
        {
            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            r += offset;
            offset = r.index();
            unsigned long limit = r.index() + a.size();
            while (r.index() < limit)
            {
                a[r.index()-offset] += *r;
                ++r;
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When adding DenseVector to DenseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                Sum<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = a.size() % parts;
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                for (int i(0); i < modulo; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div+1, i*(div+1));
                    DenseVectorRange<DT2_> range_2(b, div+1, i*(div+1));
                    TwoArgWrapper<Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div, modulo+(i*div));
                    DenseVectorRange<DT2_> range_2(b, div, modulo+(i*div));
                    TwoArgWrapper<Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(0); i < parts; ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When adding DenseVector to SparseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long modulo = b.used_elements() % parts;
            unsigned long div = b.used_elements() / parts;
            if (div == 0) 
            {
                Sum<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
                unsigned long offset;
                for (int i(0); i < modulo; ++i) 
                {
                    offset = r.index();
                    r += div;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCSum<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>,
                        const unsigned long > mywrapper(range, b, (i*(div+1)));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCSum<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>,
                        const unsigned long > mywrapper(range, b, modulo + (i*div));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i = 0; i < parts;  ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to BandedMatrix (MultiCore):");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[2*a.rows()-1];
            int taskcount(0);
            typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                {
                    TwoArgWrapper< Sum<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(*l, *r);
                    pt[taskcount] = p->dispatch(mywrapper);
                    ++taskcount; 
                }
                else
                    a.band(r.index()) = r->copy();
            }
            for (unsigned long j = 0; j < taskcount; ++j)
            {
                pt[j]->wait_on();
                delete pt[j];
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        { 
            CONTEXT("When adding DenseMatrix to DenseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
        }

        template <typename DT1_>
        static DenseMatrix<DT1_> & value(const DT1_  a, DenseMatrix<DT1_> & b)
        { 
            CONTEXT("When adding scalar to DenseMatrix (MultiCore):");

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[b.rows()];
            for (unsigned long i = 0 ; i < b.rows() ; ++i)
            {
                TwoArgWrapper< Sum<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT1_> > mywrapper(a, b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < b.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return b;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        { 
            CONTEXT("When adding DenseMatrix to SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        { 
            CONTEXT("When adding SparseMatrix to SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Sum<typename Tag_::DelegateTo>, SparseVector<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to DenseMatrix (MutiCore):");

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[2];

            ThreeArgWrapper< MCSum<typename Tag_::DelegateTo>,  DenseMatrix<DT1_>,
                const BandedMatrix<DT2_>, const bool> mywrapper1 (a, b, true);
            pt[0] = p->dispatch(mywrapper1);

            ThreeArgWrapper< MCSum<typename Tag_::DelegateTo>,  DenseMatrix<DT1_>,
                const BandedMatrix<DT2_>, const bool> mywrapper2 (a, b, false);
            pt[1] = p->dispatch(mywrapper2);

            pt[0]->wait_on();
            delete pt[0];
            pt[1]->wait_on();
            delete pt[1];
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b, const bool upper)
        {
            CONTEXT("When partial adding BandedMatrix to DenseMatrix:");

            unsigned long size(b.size());

            if (upper)
            {

                for (typename BandedMatrix<DT2_>::ConstVectorIterator r(b.band_at(size-1)), r_end(b.end_bands()) ;
                        r != r_end ; ++r)
                {
                    if (! r.exists())
                        continue;

                    unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                    unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));

                    typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                    for ( ; c != c_end ; ++c)
                    {
                        if (row_index >= size)
                            break;

                        if (col_index >= size)
                            break;

                        a[row_index][col_index] += *c;
                        ++row_index;
                        ++col_index;
                    }
                }
                return a;
            } else {
                for (typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.band_at(size-1)) ;
                        r != r_end ; ++r)
                {
                    if (! r.exists())
                        continue;

                    unsigned long size(b.size());
                    unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                    unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));

                    typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                    if (r.index() < size - 1)
                    {
                            c += ((size-1) - r.index());
                    }

                    for ( ; c != c_end ; ++c)
                    {
                        a[row_index][col_index] += *c;
                        ++row_index;
                        ++col_index;
                    }
                }
                return a;
            }
        }
    };
}
#endif
