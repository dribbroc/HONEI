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

#ifndef LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH
#define LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH 1

#include <libutil/memory_backend.hh>

namespace honei
{
    /**
     * CellBackend is the MemoryBackend for Cell-based calculations.
     *
     * \ingroup grpbackends
     * \ingroup grpgpubackend
     */
    class CellBackend :
        public MemoryBackend
    {
        public:
            /// Chunk in SPE local store memory.
            struct Chunk;

        private:
            struct Implementation;

            /// Our private implementation data.
            Implementation * _imp;

            /// Constructor.
            CellBackend();

            /// Destructor.
            ~CellBackend();

        public:
            friend class Chunk;

            /// Return the only instance of CellBackend.
            static CellBackend * instance();

            /**
             * Look up the memory chunk that is associated with a given memory id.
             *
             * \param id Memory id whose remote memory chunk shall be looked up.
             * \param device Id of the device on which the chunk exists.
             */
            const CellBackend::Chunk * get_chunk(const MemoryId id, const DeviceId id);

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
            virtual void upload(const MemoryId, const DeviceId device, void * address, const std::ptrdiff_t size);

            /**
             * Download a memory chunk from remote memory to local memory at a custom address
             *
             * \param id Memory id that uniquely identifies a remove memory chunk.
             * \param id Id of the device where to copy from.
             * \param address Local memory address whence to copy.
             * \param size Size of the memory block that will be copied.
             */
            virtual void download(const MemoryId, const DeviceId device, void * address, const std::ptrdiff_t size);

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             * \param device Id of the device on which memory shall be freed.
             */
            virtual void free(const MemoryId id, const DeviceId device);

            /**
             * Swap all memory information of two memory ids.
             *
             * \warning Both ids need to describe memory chunks of identical size.
             *
             * \param left One of the memory ids that shall be swapped.
             * \param right idem
             */
            virtual void swap(const MemoryId left, const MemoryId right);

            /// \}

            /**
             * Allocate an anonymous (i.e not associated with any memory id) memory chunk.
             *
             * \param size Size of the allocated chunk.
             */
            const CellBackend::Chunk * alloc(const DeviceId device, const unsigned int size);

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
             * Free an anonymous memory chunk.
             *
             * \param chunk Chunk that shall be freed.
             */
            void free(const CellBackend::Chunk * chunk);
    };

    /**
     * \brief Chunk of Cell SPE local store memory.
     *
     * CellBackend::Chunk represents a memory chunk in local store memory.
     * It is unique per spe context and address and can only be created, copied and freed
     * by CellBackend.
     *
     * \ingroup grpgpubackend.
     */
    struct CellBackend::Chunk
    {
        public:
            /// The id that identifies our SPE.
            const DeviceId id;

            /// Our address in the SPE's local store.
            const unsigned int address;

            /// Our size.
            const unsigned int size;

            friend class CellBackend;
            friend class CellBackend::Implementation;
        private:
            /// Constructor.
            Chunk(const DeviceId i, unsigned int a, unsigned int s);

            /// Copy-constructor.
            Chunk(const Chunk & other);

            /// Destructor.
            ~Chunk();
    };

    inline std::ostream & operator<< (std::ostream & lhs, const CellBackend::Chunk & c)
    {
        lhs << "[" << c.address << "-" << c.address + c.size << ") @" << c.id << std::endl;

        return lhs;
    }

}

#endif
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH
#define LIBUTIL_GUARD_MEMORY_BACKEND_CELL_HH 1

#include <libutil/memory_backend.hh>

namespace honei
{
    /**
     * CellBackend is the MemoryBackend for Cell-based calculations.
     *
     * \ingroup grpbackends
     * \ingroup grpcellbackend
     */
    class CellBackend :
        public MemoryBackend
    {
        public:
            /// Chunk in Cell memory.
            struct Chunk;

        private:
            struct Data;

            /// Our private implementation data.
            CellBackend::Data * _data;

            /// Constructor.
            CellBackend();

            /// Destructor.
            ~CellBackend();

            /**
             * Bind a free framebuffer object to a given Chunk.
             *
             * \param chunk Chunk that shall be associated with a framebuffer object.
             */
            void _bind_framebuffer_object(Chunk * chunk);

            /**
             * Release a framebuffer object that is associated with a given Chunk.
             *
             * \param chunk Chunk whose binding to a framebuffer object shall be released.
             */
            void _release_framebuffer_object(Chunk * chunk);

        public:
            friend Chunk;

            /// Return the only instance of CellBackend.
            static CellBackend * instance();

            /**
             * Look up the memory chunk that is associated with a given memory id.
             *
             * \param id Memory id whose remote memory chunk shall be looked up.
             */
            const CellBackend::Chunk & get_chunk(const MemoryId id);

            /// \name MemoryBackend interface
            /// \{

            /// Return the only instance of CellBackend, downcasted to MemoryBackend.
            static MemoryBackend * backend_instance();

            /**
             * Upload a memory chunk from local memory to remote memory.
             *
             * \param id Associated memory id. Valid value will be filled in if left zero.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory chunk that will be copied.
             */
            virtual void upload(const MemoryId, const void * address, const std::ptrdiff_t size);

            /**
             * Download a memory chunk from remote memory to local memory at a custom address
             *
             * \param id Memory id that uniquely identifies a remove memory chunk.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory block that will be copied. \todo eliminate?
             */
            virtual void download(const MemoryId, const void * address, const std::ptrdiff_t size);

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             */
            virtual void free(const MemoryId id);

            /// \}

            /**
             * Allocate an anonymous (i.e not associated with any memory id) memory chunk.
             *
             * \param size Size of the allocated chunk.
             */
            const CellBackend::Chunk * alloc(unsigned long size);

            /**
             * Free an anonymous memory chunk.
             *
             * \param chunk Chunk that shall be freed.
             */
            void free(const CellBackend::Chunk * chunk);
    };

    /**
     * \brief Chunk of Cell device memory.
     *
     * CellBackend::Chunk represents a memory chunk in Cell memory.
     * It is unique per effective address and can only be created, copied and freed
     * by CellBackend.
     *
     * \ingroup grpcellbackend.
     */
    struct CellBackend::Chunk
    {
        public:
            /// Our unique identifier of a Cell memory block, the effective
            /// address.
            unsigned long long ea;

            /// Our size.
            const unsigned long size;

            friend class CellBackend;

        private:
            /// Constructor.
            Chunk(unsigned long our_size);

            /// Copy-constructor.
            Chunk(const Chunk & other);

            /// Destructor.
            ~Chunk();
    };
}

#endif
