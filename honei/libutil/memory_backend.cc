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

#include <honei/libutil/memory_backend.hh>
#include <honei/libutil/stringify.hh>

using namespace honei;

MemoryBackendError::MemoryBackendError(const tags::TagValue tag, const DeviceId device, const std::string & msg) :
    Exception(msg + " in '" + stringify(tag) + "'-backend, device '" + stringify(device) + "'")
{
}

OutOfMemoryError::OutOfMemoryError(const tags::TagValue tag, const DeviceId device) :
    MemoryBackendError(tag, device, "Out of Memory")
{
}

MisalignmentError::MisalignmentError(const void * address, const unsigned alignment, const tags::TagValue tag,
        const DeviceId device) :
    MemoryBackendError(tag, device, "Misaligned address '" + stringify(address) + "', expected '" +
            stringify(alignment) + "'-byte-alignment")
{
}
