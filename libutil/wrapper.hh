/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Joachim Messer <joachim.messer@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_WRAPPER_HH
#define LIBLA_GUARD_WRAPPER_HH 1

#include <libutil/lock.hh>
#include <libutil/mutex.hh>

namespace honei
{
    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_>
        class TwoArgWrapper
        {
            public:

                typedef void result_type;

                TwoArgWrapper()
                {
                }

                void operator()(DT1_ * result, DT2_ & a, const DT3_ & b)
                {
                    *result = Tag_::value(a, b);
                }
                
                void operator()(DT1_ * result, DT2_ & a, const DT3_ & b, Mutex * mutex)
                {
                    DT1_ temp = Tag_value(a,b);
                    Lock l(*mutex);
                    *result += temp;
                }
                
        };
}
#endif


