/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBLA_GUARD_SCALED_SUM_MC_HH
#define LIBLA_GUARD_SCALED_SUM_MC_HH 1

#include <honei/libutil/partitioner.hh>
#include <honei/libutil/pool_task.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libutil/wrapper.hh>

#include <list>
#include <tr1/functional>

namespace honei
{
    // Forward declaration.
    template <typename Tag_> struct ScaledSum;

    template <typename Tag_> struct MCScaledSum
    {
        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a,
                const DenseVectorContinuousBase<DT2_> & b, const DT2_ c)
        {
            CONTEXT("In ScaledSum(DVCB, DVCB, Scalar) (MultiCore): ");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(a.size(), b.size());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scaled-sum[DVCB,DVCB,DT]::min-part-size", 1024));

            if (a.size() < 2 * min_part_size)
            {
                ScaledSum<typename Tag_::DelegateTo>::value(a, b, c);
            }
            else
            {
                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>(
                        Configuration::instance()->get_value("mc::scaled-sum[DVCB,DVCB,DT]::max-count",
                            Configuration::instance()->get_value("mc::number-of-threads", 2)),
                        min_part_size, 16, a.size(), PartitionList::Filler(partitions));

                for (PartitionList::ConstIterator p(partitions.begin()),
                        p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    DenseVectorRange<DT1_> a_range(a.range(p->size, p->start));
                    DenseVectorRange<DT2_> b_range(b.range(p->size, p->start));

                    ThreeArgWrapper< ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                        const DenseVectorRange<DT2_>, const DT2_> wrapper(a_range, b_range, c);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a,
                const DenseVectorContinuousBase<DT2_> & b, const DenseVectorContinuousBase<DT2_> & c)
        {
            CONTEXT("In ScaledSum(DVCB, DVCB, DVCB) (MultiCore): ");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(a.size(), b.size());
            if (a.size() != c.size())
                throw VectorSizeDoesNotMatch(a.size(), c.size());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scaled-sum[DVCB,DVCB,DVCB]::min-part-size", 1024));

            if (a.size() < 2 * min_part_size)
            {
                ScaledSum<typename Tag_::DelegateTo>::value(a, b, c);
            }
            else
            {
                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>(
                        Configuration::instance()->get_value("mc::scaled-sum[DVCB,DVCB,DVCB]::max-count",
                            Configuration::instance()->get_value("mc::number-of-threads", 2)),
                        min_part_size, 16, a.size(), PartitionList::Filler(partitions));

                for (PartitionList::ConstIterator p(partitions.begin()),
                        p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    DenseVectorRange<DT1_> a_range(a.range(p->size, p->start));
                    DenseVectorRange<DT2_> b_range(b.range(p->size, p->start));
                    DenseVectorRange<DT2_> c_range(c.range(p->size, p->start));

                    ThreeArgWrapper< ScaledSum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>,
                        const DenseVectorRange<DT2_>, const DenseVectorRange<DT2_> >
                            wrapper(a_range, b_range, c_range);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }

            return a;
        }
    };

}

#endif
