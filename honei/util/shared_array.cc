/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/util/shared_array-impl.hh>

namespace honei
{
    SharedArrayError::SharedArrayError(const std::string & message) throw () :
        Exception(message)
    {
    }

    template class SharedArray<float>;

    template class SharedArray<double>;

    template class SharedArray<int>;

    template class SharedArray<unsigned int>;

    template class SharedArray<long>;

    template class SharedArray<unsigned long>;

    template class SharedArray<bool>;

#ifdef HONEI_GMP
    template class SharedArray<mpf_class>;
#endif
}
