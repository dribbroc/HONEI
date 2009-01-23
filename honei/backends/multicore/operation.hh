/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Sven Mallach <sven.mallach@cs.uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MULTICORE_GUARD_OPERATION_HH
#define MULTICORE_GUARD_OPERATION_HH 1

#include <honei/backends/multicore/dispatch_policy.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/operation_wrapper.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/tags.hh>

#include <tr1/functional>

namespace honei
{
    namespace mc
    {
        template <typename DelegateOperationType_> struct Operation
        {
            template <typename DT1_, typename DT2_>
            static void op(DenseVectorBase<DT1_> & x, const DT2_ & a, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, a);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());

                        OperationWrapper<DelegateOperationType_, DenseVectorBase<DT1_>,
                        DenseVectorBase<DT1_>, DenseVectorBase<DT2_>, DT2_> wrapper(x_slice);

                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_slice, a)));
                    }

                    DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());

                    DelegateOperationType_::value(x_slice, a);

                    tickets.wait();
                }

                return x;
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorContinuousBase<DT1_> & x, const DT2_ & a, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, a);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));

                        OperationWrapper<DelegateOperationType_, DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT1_>, DT2_ > wrapper(x_range);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_range, a)));
                    }

                    DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));

                    DelegateOperationType_::value(x_range, a);

                    tickets.wait();
                }

                return x;
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                        DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() + (p->start * y.stepsize()), y.stepsize());

                        OperationWrapper<DelegateOperationType_, DenseVectorBase<DT1_>, DenseVectorBase<DT1_>, DenseVectorBase<DT2_> > wrapper(x_slice);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_slice, y_slice)));
                    }

                    DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                    DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() + (p->start * y.stepsize()), y.stepsize());

                    DelegateOperationType_::value(x_slice, y_slice);

                    tickets.wait();
                }
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                        DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));

                        OperationWrapper<DelegateOperationType_, DenseVectorContinuousBase<DT1_>,
                        DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT2_> > wrapper(x_range);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_range, y_range)));
                    }

                    DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                    DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));

                    DelegateOperationType_::value(x_range, y_range);

                    tickets.wait();
                }
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y, DT2_ b, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y, b);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                        DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() + (p->start * y.stepsize()), y.stepsize());

                        OperationWrapper<DelegateOperationType_, DenseVectorBase<DT1_>,
                        DenseVectorBase<DT1_>, DenseVectorBase<DT2_>, DT2_> wrapper(x_slice);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_slice, y_slice, b)));
                    }

                    DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                    DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() + (p->start * y.stepsize()), y.stepsize());

                    DelegateOperationType_::value(x_slice, y_slice, b);

                    tickets.wait();
                }
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y, DT2_ b,
                    unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y, b);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                        DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));

                        OperationWrapper<DelegateOperationType_, DenseVectorContinuousBase<DT1_>,
                        DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT2_>, DT2_> wrapper(x_range);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_range, y_range, b)));
                    }

                        DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                        DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));

                    DelegateOperationType_::value(x_range, y_range, b);

                    tickets.wait();
                }
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y, const DenseVectorBase<DT2_> & z,
                    unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() != z.size())
                    throw VectorSizeDoesNotMatch(z.size(), x.size());


                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y, z);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                        DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() +
                        (p->start * y.stepsize()), y.stepsize());
                        DenseVectorSlice<DT1_> z_slice(z, p->size, z.offset() +
                        (p->start * z.stepsize()), z.stepsize());

                        OperationWrapper<DelegateOperationType_, DenseVectorBase<DT1_>,
                            DenseVectorBase<DT1_>, DenseVectorBase<DT2_>,
                            DenseVectorBase<DT2_> > wrapper(x_slice);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper,
                        x_slice, y_slice, z_slice)));
                    }

                    DenseVectorSlice<DT1_> x_slice(x, p->size, x.offset() + (p->start * x.stepsize()), x.stepsize());
                    DenseVectorSlice<DT1_> y_slice(y, p->size, y.offset() +
                        (p->start * y.stepsize()), y.stepsize());
                    DenseVectorSlice<DT1_> z_slice(z, p->size, z.offset() +
                        (p->start * z.stepsize()), z.stepsize());

                    DelegateOperationType_::value(x_slice, y_slice, z_slice);

                    tickets.wait();
                }
            }

            template <typename DT1_, typename DT2_>
            static void op(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y,
                    const DenseVectorContinuousBase<DT2_> & z, unsigned long min_part_size, unsigned long max_count)
            {
                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                if (x.size() != z.size())
                    throw VectorSizeDoesNotMatch(z.size(), x.size());

                if (x.size() < 2 * min_part_size)
                {
                    DelegateOperationType_::value(x, y, z);
                }
                else
                {
                    PartitionList partitions;
                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                    TicketVector tickets;

                    PartitionList::ConstIterator p(partitions.begin());
                    for (PartitionList::ConstIterator p_last(partitions.last()) ; p != p_last ; ++p)
                    {
                        DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                        DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));
                        DenseVectorRange<DT2_> z_range(z.range(p->size, p->start));

                        OperationWrapper<DelegateOperationType_, DenseVectorContinuousBase<DT1_>,
                            DenseVectorContinuousBase<DT1_>, DenseVectorContinuousBase<DT2_>,
                            DenseVectorContinuousBase<DT2_> > wrapper(x_range);
                        tickets.push_back(mc::ThreadPool::instance()->enqueue(std::tr1::bind(wrapper, x_range, y_range, z_range)));
                    }

                    DenseVectorRange<DT1_> x_range(x.range(p->size, p->start));
                    DenseVectorRange<DT2_> y_range(y.range(p->size, p->start));
                    DenseVectorRange<DT2_> z_range(z.range(p->size, p->start));

                    DelegateOperationType_::value(x_range, y_range, z_range);

                    tickets.wait();
                }
            }
        };
    }
}
#endif
