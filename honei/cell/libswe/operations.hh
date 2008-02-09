/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CELL_GUARD_LIBSWE_OPERATIONS_HH
#define CELL_GUARD_LIBSWE_OPERATIONS_HH 1

#include <honei/cell/libutil/operation_framework.hh>

namespace honei
{
    namespace cell
    {
        namespace operations
        {
            extern Operation<1, float, rtm_dma> flow_processing_x_float;

            extern Operation<1, float, rtm_dma> flow_processing_y_float;
        }
    }
}

#endif
