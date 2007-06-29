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

#ifndef LIBUTIL_GUARD_MEMORY_MANAGER_HH
#define LIBUTIL_GUARD_MEMORY_MANAGER_HH 1

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/memory_backend.hh>
#include <libutil/tags.hh>
#include <libutil/stringify.hh>

#include <map>
#include <set>

namespace pg512
{
    /**
     * MemoryId uniquly identifies remote memory chunks across all backends.
     *
     * \ingroup grpmemorymanager
     */
    typedef unsigned long long MemoryId;

    /**
     * MemoryAddressNotKnown is thrown by MemoryManager during address lookup.
     *
     * \ingroup grpexceptions
     * \ingroup grpmemorymanager
     */
    struct MemoryAddressNotKnown :
        public Exception
    {
        /**
         * Constructor.
         *
         * \param address Unknown local memory address.
         */
        MemoryAddressNotKnown(const void * address);
    };

    /**
     * MemoryIdNotKnown is thrown by MemoryManager when it encounters an unused
     * MemoryId in transfers and free.
     *
     * \ingroup grpexceptions
     * \ingroup grpmemorymanager
     */
    struct MemoryIdNotKnown :
        public Exception
    {
        /**
         * Constructor.
         *
         * \param id Unknown memory id.
         */
        MemoryIdNotKnown(const MemoryId id);
    };

    /**
     * MemoryChunkSizeInvalid is thrown by memory backends whenever the requested
     * chunk size is illegal for the particular backend.
     *
     * \ingroup grpexceptions
     * \ingroup grpmemorybackends
     */
    struct MemoryChunkSizeInvalid :
        public Exception
    {
        /**
         * Constructor.
         *
         * \param size Requested size of the memory chunk.
         * \param msg Error message that explains why size is illegal for the backend.
         */
        MemoryChunkSizeInvalid(unsigned long size, const std::string & msg);
    };

    /**
     * MemoryInfo holds all information that is associated with a MemoryId.
     *
     * \ingroup grpmemorymanager
     */
    class MemoryInfo
    {
        private:
            /// Our address in local memory.
            void * _address;

            /// Constructor.
            MemoryInfo(void * a, std::ptrdiff_t s, tags::TagValue l) :
                _address(a),
                location(l),
                size(s),
                ttime(clock())
            {
            }

        public:
            friend class MemoryManager;
            /// Our remote location's namespace.
            const tags::TagValue location;

            /// Our size of local memory.
            const std::ptrdiff_t size;

            /// Our time of last transfer.
            const clock_t ttime;

            /// Returns our address in local memory.
            inline const void * address() const
            {
                return _address;
            }

            /// Copy-constructor.
            MemoryInfo(const MemoryInfo & other) :
                _address(other._address),
                location(other.location),
                size(other.size),
                ttime(other.ttime)
            {
            }
    };

    /**
     * MemoryManager is used to transfer memory chunks from local memory to remote memory and vice versa.
     * It does not implement any particular ways to transfer memory. Rather it uses descendants of
     * MemoryBackend to accomplish these tasks.
     *
     * \ingroup grpmemorymanager
     */
    class MemoryManager
    {
        public:
            /// \name Convenience typedef for backend handling
            /// \{
            typedef MemoryBackend * (*InstanceFunction)();
            typedef std::map<tags::TagValue, InstanceFunction> BackendMap;
            /// \}

        private:
            /// \name Convenience types for memory handling
            /// \{
            typedef std::set<MemoryId> IdSet;
            typedef std::map<MemoryId, MemoryInfo *> InfoMap;
            typedef std::map<MemoryId, void *> AddressMap;
            /// \}

            /// Our set of memory ids. Unique by timestamp.
            IdSet _ids;

            /// Our map of memory ids to memory information. Unique per id.
            InfoMap _info_map;

            /// Our map of memory ids to local memory addresses. Unique per memory address.
            AddressMap _address_map;

            /// Our map of tag values to backend instance functions. Unique per tag value.
            BackendMap _backend_map;

            /// \name Basic operations
            /// \{
            /// Constructor.
            MemoryManager();

            /// Destructor.
            ~MemoryManager();

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            MemoryManager(const MemoryManager &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            const MemoryManager & operator= (const MemoryManager &);
            /// \}

            /// Return true if a memory id is in use.
            inline bool _id_used(const MemoryId id) const
            {
                return _ids.end() != _ids.find(id);
            }

            /// Return the next free memory id.
            inline MemoryId _get_free_id()
            {
                static MemoryId last_id(0);

                do
                {
                    ++last_id;
                }
                while (_id_used(last_id) || (last_id == 0));

                _ids.insert(last_id);
                return last_id;
            }

            /// Return the backend for a given destination.
            inline MemoryBackend * _get_backend(const tags::TagValue destination)
            {
                BackendMap::const_iterator b(_backend_map.find(destination));
                ASSERT(_backend_map.end() != b, "No backend found for tag value '" +
                        stringify(destination) + "'");

                return b->second();
            }

        public:
            friend class MemoryBackendRegistrator;

            /// Return the only instance of MemoryManager.
            static MemoryManager * instance();

            /**
             * Upload a memory chunk from local memory to remote memory.
             *
             * \param id Associated memory id. Valid value will be filled in if left zero.
             * \param destination Tag value of the remote memory whence to copy to.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory chunk that will be copied.
             */
            void upload(MemoryId & id, const tags::TagValue location, void * address, const std::ptrdiff_t size);

            /**
             * Download a memory chunk from remote memory to local memory.
             *
             * \param id Memory id that uniquely identifies a remote memory chunk.
             */
            void download(const MemoryId id);

            /**
             * Download a memory chunk from remote memory to local memory at a custom address
             *
             * \param id Memory id that uniquely identifies a remove memory chunk.
             * \param address Local memory address whence to copy from.
             * \param size Size of the memory block that will be copied.
             */
            void download(const MemoryId id, void * address, const std::ptrdiff_t size);

            /**
             * Free an existing memory id and its associated remote memory.
             * Local memory will not be tempered with.
             *
             * \param id Memory id that shall be freed.
             */
            void free(const MemoryId id);

            /**
             * Look up the memory id of a given local memory address.
             *
             * \param address Local memory address that shall be looked up.
             */
            const MemoryId get_id_by_address(const void * address) const;

            /**
             * Swap all memory information of two memory ids.
             *
             * \warning Both ids need to describe memory chunks of identical size.
             *
             * \param left One of the memory ids that shall be swapped.
             * \param right idem
             */
            void swap(const MemoryId left, const MemoryId right);
    };

    /**
     * MemoryBackendRegistrator registers a descendant of MemoryBackend with MemoryManager.
     *
     * \ingroup grpmemorymanager
     */
    struct MemoryBackendRegistrator
    {
        /**
         * Constructor.
         *
         * \param v Tagvalue that the backend is associated with.
         * \param f Singleton-instance function of the backend.
         */
        MemoryBackendRegistrator(const tags::TagValue v, MemoryManager::InstanceFunction f)
        {
            MemoryManager::instance()->_backend_map.insert(std::make_pair(v, f));
        }
    };
}

#endif
