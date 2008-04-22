/* vim: set sw=4 sts=4 et nofoldenable : */
/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Joachim Messer <joachim.messer@uni-dortmund.de>
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

#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>

#include <list>

namespace honei
{
    template <> struct Implementation<PartitionList>
    {
        std::list<Partition> partitions;
    };
}

using namespace honei;

PartitionList::PartitionList() :
    PrivateImplementationPattern<PartitionList, Shared>(new Implementation<PartitionList>)
{
}

PartitionList::~PartitionList()
{
}

PartitionList::ConstIterator
PartitionList::begin() const
{
    return ConstIterator(_imp->partitions.begin());
}

PartitionList::ConstIterator
PartitionList::last() const
{
    return ConstIterator(--_imp->partitions.end());
}

PartitionList::ConstIterator
PartitionList::end() const
{
    return ConstIterator(_imp->partitions.end());
}

unsigned
PartitionList::size() const
{
    return _imp->partitions.size();
}

void
PartitionList::Filler::operator() (unsigned long start, unsigned long size)
{
    _partition_list._imp->partitions.push_back(Partition(start, size));
}

Partitioner<tags::Cell>::Partitioner(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size,
    std::tr1::function<void(unsigned long, unsigned long)> dispatch)
{
    CONTEXT("When partitioning problem of size '" + stringify(overall_size) + "' (Cell):");

    unsigned part_size(0);
    unsigned rest(0);
    unsigned count(0);

    if ((best_part_size > overall_size))
    {
        count = 1;
        part_size = overall_size;
    }
    else if (max_count * best_part_size >= overall_size)
    {
        count = overall_size / best_part_size;
        rest = overall_size % best_part_size;
        rest = rest - rest % 16;
        part_size = best_part_size - best_part_size % 16;
    }
    else
    {
        part_size = overall_size / max_count;
        part_size = part_size - part_size % 16;
        rest = overall_size % part_size;
        rest = rest - rest % 16;
        count = overall_size / part_size;
    }

    if (part_size > 0)
    {
        dispatch(0, part_size + rest);
    }

    unsigned long start(part_size + rest);

    for (unsigned i(1); i < count; ++i)
    {
        if (part_size > 0) 
        {
            dispatch(start, part_size);
        }

        start += part_size;
    }

    if (overall_size > start)
    {
        dispatch(start, overall_size - start);
    }
}

Partitioner<tags::CPU::MultiCore>::Partitioner(unsigned long max_count, unsigned long best_part_size,
        unsigned long quantization, unsigned long overall_size,
        std::tr1::function<void(unsigned long, unsigned long)> dispatch)
{
    CONTEXT("When partitioning problem of size '" + stringify(overall_size) + "' (MC):");

    unsigned part_size(0);
    unsigned count(0);
    unsigned modulo(0);

    if (best_part_size < quantization)
    {
        best_part_size = quantization;
    }
    if (overall_size < (best_part_size << 1))
    {
        dispatch(0, overall_size);
    }
    else
    {
        count = overall_size / best_part_size;

        if (count > max_count)
        {
            count = max_count;
        }

            part_size = overall_size / count;
            part_size = part_size - part_size % quantization;
            modulo = (overall_size - count * part_size) / quantization;

        unsigned long start(0);

        for (unsigned i(0); i < modulo; ++i)
        {
            dispatch(start, part_size + quantization);
            start += part_size + quantization;
        }

        for (unsigned i(modulo); i < count - 1; ++i)
        {
            dispatch(start, part_size);
            start += part_size;
        }
        dispatch(start, overall_size - start);
    }
}

