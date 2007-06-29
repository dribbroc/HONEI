/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
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

#ifndef LIBUTIL_GUARD_MEMORY_BACKEND_GPU_HH
#define LIBUTIL_GUARD_MEMORY_BACKEND_GPU_HH 1

#include <libutil/memory_backend.hh>

#include <GL/glew.h>

namespace pg512
{
    /**
     * GPUBackend is the MemoryBackend for GPU-based calculations.
     *
     * \ingroup grpbackends
     * \ingroup grpgpubackend
     */
    class GPUBackend :
        public MemoryBackend
    {
        public:
            /// Chunk in GPU memory.
            struct Chunk;

        private:
            struct Data;

            /// Our private implementation data.
            GPUBackend::Data * _data;

            /// Constructor.
            GPUBackend();

            /// Destructor.
            ~GPUBackend();

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
            friend class Chunk;

            /// Return the only instance of GPUBackend.
            static GPUBackend * instance();

            /**
             * Look up the memory chunk that is associated with a given memory id.
             *
             * \param id Memory id whose remote memory chunk shall be looked up.
             */
            const GPUBackend::Chunk & get_chunk(const MemoryId id);

            /// \name MemoryBackend interface
            /// \{

            /// Return the only instance of GPUBackend, downcasted to MemoryBackend.
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
            virtual void download(const MemoryId, void * address, const std::ptrdiff_t size);

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             */
            virtual void free(const MemoryId id);

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
            const GPUBackend::Chunk * alloc(unsigned long size);

            /**
             * Free an anonymous memory chunk.
             *
             * \param chunk Chunk that shall be freed.
             */
            void free(const GPUBackend::Chunk * chunk);
    };

    /**
     * \brief Chunk of GPU device memory.
     *
     * GPUBackend::Chunk represents a memory chunk in GPU memory.
     * It is unique per texture id and can only be created, copied and freed
     * by GPUBackend.
     *
     * \ingroup grpgpubackend.
     */
    struct GPUBackend::Chunk
    {
        public:
            /// Our unique identifier of a GPU memory block.
            GLuint id;

            /// Our size.
            const unsigned long size;

            /// Our width (internal representation).
            const GLuint width;

            /// Our associated framebuffer object.
            GLuint fbo;

            /// Our associated texture target.
            const static GLenum texture_target = GL_TEXTURE_RECTANGLE_ARB;

            /// \todo
            const static GLenum internal_format = GL_FLOAT_R32_NV;

            /// \todo
            const static GLenum format = GL_LUMINANCE;

            /// Bind free framebuffer object to us.
            inline void bind_framebuffer_object()
            {
                GPUBackend::instance()->_bind_framebuffer_object(this);
            }

            /// Release our framebuffer object;
            inline void release_framebuffer_object()
            {
                GPUBackend::instance()->_release_framebuffer_object(this);
            }

            friend class GPUBackend;

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
