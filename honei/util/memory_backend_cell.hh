/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Dirk Ribbrock <d_ribbrock@web.de>
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

#ifndef LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH
#define LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/memory_backend.hh>

#include <vector>

namespace honei
{
    /**
     * CellBackend is the MemoryBackend for Cell-based calculations.
     *
     * \ingroup grpbackends
     * \ingroup grpgpubackend
     */
    class CellBackend :
        public MemoryBackend,
        public InstantiationPolicy<CellBackend, Singleton>
    {
       private:
            struct Implementation;
            struct AccessCounter;

            typedef std::vector<AccessCounter> SPECounters;

            /// Our private implementation data.
            Implementation * _imp;

            /// Our SPE Counter Vector
            SPECounters _spe_counters;

            /// Constructor.
            CellBackend();

        public:
            friend class InstantiationPolicy<CellBackend, Singleton>;

            /// \name MemoryBackend interface
            /// \{

            /// Return the only instance of CellBackend, downcasted to MemoryBackend.
            static MemoryBackend * backend_instance();

            /**
             * Upload a memory chunk from local memory to remote memory.
             *
             * \param id Associated memory id. Valid value will be filled in if left zero.
             * \param device Id of the device whence to copy.
             * \param address Local memory address where to copy from.
             * \param size Size of the memory chunk that will be copied.
             */
            virtual void upload(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size);

            /**
             * Download a memory chunk from remote memory to local memory at a custom address
             *
             * \param id Memory id that uniquely identifies a remove memory chunk.
             * \param device Id of the device where to copy from.
             * \param address Local memory address whence to copy.
             * \param size Size of the memory block that will be copied.
             */
            virtual void download(const MemoryId id, const DeviceId device, void * address, const std::ptrdiff_t size);

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             * \param device Id of the device on which memory shall be freed.
             */
            virtual void free(const MemoryId id, const DeviceId device);

            /**
             * Swap all memory information of two memory ids on all devices of a
             * backend.
             *
             * \warning Both ids need to describe memory chunks of identical size.
             *
             * \param left One of the memory ids that shall be swapped.
             * \param right idem
             */
            virtual void swap(const MemoryId left, const MemoryId right);

            /// \}

            /**
             * Swap all memory information of two memory ids for a specific
             * SPE only.
             *
             * \warning Both ids need to describe memory chunks of identical size.
             *
             * \param left One of the memory ids that shall be swapped.
             * \param right idem
             * \param device Id of the device.
             */
            virtual void swap(const MemoryId left, const MemoryId right, const DeviceId device);

            /**
             * Mark device as holding data for reading.
             *
            **/
            void read_enable(DeviceId device_id);

            /**
             * Unmark device as holding data for reading.
             *
            **/
            void read_disable(DeviceId device_id);

            /**
             * Mark device as holding data for writing.
             *
            **/
            void write_enable(DeviceId device_id);

            /**
             * Unmark device as holding data for writing.
             *
            **/
            void write_disable(DeviceId device_id);

   };

}

#endif

