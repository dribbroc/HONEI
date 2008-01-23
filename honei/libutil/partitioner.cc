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

#include <honei/libutil/exception.hh>
#include <honei/libutil/partitioner.hh>

#include <tr1/functional>

using namespace honei;

Partitioner::Partitioner(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size,
        std::tr1::function<void(unsigned long, unsigned long)> dispatch)
{
    CONTEXT("When partitioning problem of size " + overall_size);
    std::list<Parts> partition(Partitioner::partition(max_count, best_part_size, overall_size));
    for (std::list<Parts>::iterator i(partition.begin()), i_end(partition.end()) ;
            i != i_end ; ++i)
    {
        dispatch(i->start, i->size);
    }
}

std::list<Parts> Partitioner::partition(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size)
{
    CONTEXT("When partitioning problem of size " + overall_size);

    std::list<Parts> result;
    unsigned part_size(0);
    unsigned count(0);
    if (best_part_size >= overall_size)
    {
        part_size = (overall_size - overall_size % 32) / 2;
        if (part_size > 0)
        {
            result.push_back(Parts(0, part_size));
            result.push_back(Parts(part_size, part_size));
        }
        if (overall_size > 2 * part_size)
        {
            result.push_back(Parts(2 * part_size, overall_size - 2 * part_size));
        }

    }
    else
    {
        if (2 * best_part_size >= overall_size)
        {
            part_size = best_part_size;
            count = 1;
        }
        else
        {
            count = overall_size/best_part_size;
            if (count > max_count)
            {
                count = max_count;
            }
            part_size = overall_size/count;
            part_size = part_size - part_size % 16;
        }

        unsigned long start(0);
        for (unsigned i(0); i < count; ++i)
        {
            if (part_size > 0) 
                result.push_back(Parts(start, part_size));
            start += part_size;
        }

        if (overall_size > start)
        {
            result.push_back(Parts(start, overall_size - start));
        }
    }
    return result;
}

