/* vim: set sw=4 sts=4 et nofoldenable : */
/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Joachim Messer <joachim.messer@uni-dortmund.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBUTIL_GUARD_PARTITIONER_HH
#define LIBUTIL_GUARD_PARTITIONER_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/tags.hh>
#include <honei/util/wrapped_forward_iterator.hh>

#include <tr1/functional>
#include <tr1/memory>

namespace honei
{
    struct Partition
    {
        unsigned long start;
        unsigned long size;

        Partition(unsigned long st, unsigned long si) :
            start(st),
            size(si)
        {
        }
    };

    class PartitionList :
        public PrivateImplementationPattern<PartitionList, Shared>
    {
        public:
            /// \name Basic operations
            /// \{

            /// Constructor.
            PartitionList();

            /// Destructor.
            ~PartitionList();

            /// \}

            struct Filler;

            /// \name Iteration over our elements
            /// \{

            struct ConstIteratorTag;
            typedef WrappedForwardIterator<ConstIteratorTag, Partition> ConstIterator;

            ConstIterator begin() const;

            ConstIterator last() const;

            ConstIterator end() const;

            /// \}

            /// Return or size.
            unsigned size() const;
    };

    class PartitionList::Filler
    {
        private:
            /// Our partition list.
            PartitionList _partition_list;

        public:
            /// Constructor.
            Filler(const PartitionList & partition_list) :
                _partition_list(partition_list)
            {
            }

            /// Evaluation operator, used to fill the PartitionList.
            void operator() (unsigned long start, unsigned long size);
    };

    /**
     * \name Partitioner
     * \{
     * Splits a problem of size overall_size into (at most) max_count + 1 parts, using at least best_part_size
     * elements per partition for all but the last partition. Dispatches each partition to dispatch.
     */

    template <typename Tag_> class Partitioner;

    template <> class Partitioner<tags::Cell> :
        public InstantiationPolicy<Partitioner<tags::Cell>, NonCopyable>
    {
        public:
            Partitioner(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size,
                    std::tr1::function<void(unsigned long, unsigned long)> dispatch);
    };

    template <> class Partitioner<tags::CPU::MultiCore> :
        public InstantiationPolicy<Partitioner<tags::CPU::MultiCore>, NonCopyable>
    {
        public:
            Partitioner(unsigned long max_count, unsigned long best_part_size, unsigned long quantization,
                    unsigned long overall_size, std::tr1::function<void(unsigned long, unsigned long)> dispatch);
    };
}

#endif
