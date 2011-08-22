/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#pragma once
#ifndef UTIL_GUARD_CLIP_BOARD_HH
#define UTIL_GUARD_CLIP_BOARD_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/instantiation_policy-impl.hh>

#include <vector>

namespace honei
{
    template <typename ContainerType_, long magic_number_ = 1>
    class ClipBoard :
        public InstantiationPolicy<ClipBoard<ContainerType_, magic_number_>, Singleton>
    {
        private:
            /// Constructor
            ClipBoard(){};

            /// Destructor
            ~ClipBoard(){};

            std::vector<ContainerType_> _list;

        public:
            friend class InstantiationPolicy<ClipBoard, Singleton>;

            void push_back(ContainerType_ & input)
            {
                _list.push_back(input);
            }

            void pop_back()
            {
                _list.pop_back();
            }

            ContainerType_ & at(unsigned long i)
            {
                return _list.at(i);
            }

            void clear()
            {
                _list.clear();
            }

            unsigned long size()
            {
                return _list.size();
            }

    };
}
#endif
