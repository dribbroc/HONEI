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

#ifndef LIBUTIL_GUARD_MEMORY_BACKEND_HH
#define LIBUTIL_GUARD_MEMORY_BACKEND_HH 1

namespace pg512
{
    typedef unsigned long long MemoryId;

    /**
     * DeviceId uniquely identifies one of the available remote memory devices that a backend
     * can interface. Possible are e.g. several Cell BE SPE Local Stores or several video memories.
     */
    typedef unsigned long DeviceId;

    /**
     * If there is only one device supported by a MemoryBackend, it must not accept requests for
     * other devices except default_device.
     */
    const DeviceId default_device(~0x0L);

    /**
     * MemoryBackend is the interface class for all memory backends that MemoryManager supports.
     *
     * \ingroup grpmemorymanager
     * \ingroup grpmemorybackends
     */
    class MemoryBackend
    {
        public:
            /**
             * Upload a memory chunk from local memory to remote memory.
             *
             * \param id Associated memory id.
             * \param device Id of the the device whence to copy to.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory chunk that will be copied.
             */
            virtual void upload(const MemoryId id, const DeviceId device, const void * address, const std::ptrdiff_t size) = 0;

            /**
             * Download a memory chunk from remote memory to local memory at a custom address
             *
             * \param id Memory id that uniquely identifies a remove memory chunk.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory block that will be copied. \todo eliminate?
             */
            virtual void download(const MemoryId id, void * address, const std::ptrdiff_t size) = 0;

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             */
            virtual void free(const MemoryId id) = 0;

            /**
             * Swap all memory information of two memory ids.
             *
             * \warning Both ids need to describe memory chunks of identical size.
             *
             * \param left One of the memory ids that shall be swapped.
             * \param right idem
             */
            virtual void swap(const MemoryId left, const MemoryId right) = 0;
    };
}

#endif
