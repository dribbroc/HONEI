/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/cell/cell.hh>
#include <honei/cell/libla/operations.hh>
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

using namespace honei::cell;

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void sum_dense_scalar_double(vector double * a, const unsigned size, const double scalar)
            {
                vector double scalars = spu_splats(scalar);
                for (unsigned i(0) ; i < size ; ++i)
                {
                    a[i] = spu_add(a[i], scalars);
                }
            }
        }

        namespace operations
        {
            Operation<1, double, rtm_dma> sum_dense_scalar_double = {
                &implementation::sum_dense_scalar_double
            };
        }
    }
}




