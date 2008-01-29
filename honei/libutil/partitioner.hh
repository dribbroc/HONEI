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

#ifndef LIBUTIL_GUARD_PARTITIONER_HH
#define LIBUTIL_GUARD_PARTITIONER_HH 1

#include <list>
#include <tr1/functional>

namespace honei
{
    struct Parts
    {
        unsigned long start;
        unsigned long size;

        Parts()
        {
        };

        Parts(unsigned long sstart, unsigned long ssize)
        {
            start = sstart;
            size = ssize;
        };

        void operator= (Parts part)
        {
            start = part.start;
            size = part.size;
        };
    };

    class Partitioner
    {
        public:
            Partitioner(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size,
                    std::tr1::function<void(unsigned long, unsigned long)> dispatch);

            static std::list<Parts> partition(unsigned long max_count, unsigned long best_part_size, unsigned long overall_size);
    };
}

#endif
