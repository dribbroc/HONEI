/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/exception.hh>
#include <honei/util/memory_backend.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/backends/cuda/transfer.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/backends/cuda/gpu_pool.hh>

namespace honei
{
    namespace cuda
    {
        class UploadTask
        {
            private:
                void * device;
                void * address;
                unsigned long bytes;
            public:
                UploadTask(void * a, void * d, unsigned long b) :
                    device(d),
                    address(a),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_upload(address, device, bytes);
                }
        };

        class DownloadTask
        {
            private:
                void * device;
                void * address;
                unsigned long bytes;
            public:
                DownloadTask(void * d, void * a, unsigned long b) :
                    device(d),
                    address(a),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_download(device, address, bytes);
                }
        };

        class AllocTask
        {
            private:
                void ** device;
                unsigned long bytes;
            public:
                AllocTask(void ** d, unsigned long b) :
                    device(d),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    *device = cuda_malloc(bytes);
                }
        };

        class FreeTask
        {
            private:
                void * device;
            public:
                FreeTask(void * d) :
                    device(d)
                {
                }

                void operator() ()
                {
                    cuda_free(device);
                }
        };

        class CopyTask
        {
            private:
                void * src;
                void * dest;
                unsigned long bytes;
            public:
                CopyTask(void * s, void * d, unsigned long b) :
                    src(s),
                    dest(d),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_copy(src, dest, bytes);
                }
        };

        class ConvertFloatDoubleTask
        {
            private:
                void * src;
                void * dest;
                unsigned long bytes;
            public:
                ConvertFloatDoubleTask(void * s, void * d, unsigned long b) :
                    src(s),
                    dest(d),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_convert_float_double(src, dest, bytes);
                }
        };

        class ConvertDoubleFloatTask
        {
            private:
                void * src;
                void * dest;
                unsigned long bytes;
            public:
                ConvertDoubleFloatTask(void * s, void * d, unsigned long b) :
                    src(s),
                    dest(d),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_convert_double_float(src, dest, bytes);
                }
        };

        class FillTask
        {
            private:
                void * address;
                unsigned long bytes;
            public:
                FillTask(void * a, unsigned long b) :
                    address(a),
                    bytes(b)
                {
                }

                void operator() ()
                {
                    cuda_fill_zero(address, bytes);
                }
        };
    }
    /**
     * \brief MemoryBackend<tags::GPU::CUDA>::Implementation is the private implementation class for
     * MemoryBackend<tags::GPU::CUDA>.
     *
     * \ingroup grpvector
     */
    template <> struct Implementation<MemoryBackend<tags::GPU::CUDA> >
    {
        struct Chunk
        {
            public:
                Chunk(void * a, void * d, unsigned long b, int did) :
                    address(a),
                    device(d),
                    bytes(b),
                    device_id(did)
            {
            }
                void * address;
                void * device;
                unsigned long bytes;
                int device_id;
        };

        std::multimap<void *, void *> _address_map;
        std::multimap<void *, int> _device_map;

        std::multimap<void *, Chunk> _id_map;

        void * upload(void * memid, void * address, unsigned long bytes)
        {
            std::pair<std::multimap<void *, void *>::iterator, std::multimap<void *, void *>::iterator> range(_address_map.equal_range(address));
            if (range.first == range.second)
            //address unknown
            {
                bool idle(cuda::GPUPool::instance()->idle());
                // running in slave thread
                if (! idle)
                {
                    void * device(this->alloc(memid, address, bytes));
                    cuda_upload(address, device, bytes);
                    return device;
                }
                // running in main thread -> switch to slave thread
                else
                {
                    void * device(this->alloc(memid, address, bytes));
                    cuda::UploadTask ut(address, device, bytes);
                    cuda::GPUPool::instance()->enqueue(ut, cuda_get_device()).wait();
                    return device;
                }
            }
            else
            {
                // address already known on some device
                std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
                for (std::multimap<void *, Chunk>::iterator map_i(range.first) ; map_i != range.second ; ++map_i)
                {
                    if (map_i->second.address == address && map_i->second.device_id == cuda_get_device())
                        return map_i->second.device;
                }
                //address not known on our current device, alloc and upload
                {
                    bool idle(cuda::GPUPool::instance()->idle());
                    // running in slave thread
                    if (! idle)
                    {
                        void * device(this->alloc(memid, address, bytes));
                        cuda_upload(address, device, bytes);
                        return device;
                    }
                    // running in main thread -> switch to slave thread
                    else
                    {
                        void * device(this->alloc(memid, address, bytes));
                        cuda::UploadTask ut(address, device, bytes);
                        cuda::GPUPool::instance()->enqueue(ut, cuda_get_device()).wait();
                        return device;
                    }
                }
            }
        }

        void download(void * memid, void * /*address*/, unsigned long /*bytes*/)
        {
            //std::multimap<void *, Chunk>::iterator map_i;
            std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
            if (range.first == range.second)
            {
                //address not known
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::download address not found!");
            }
            else
            {
                for (std::multimap<void *, Chunk>::iterator map_i(range.first) ; map_i != range.second ; ++map_i)
                {
                    //std::map<void *, void *>::iterator i(_address_map.find(address));
                    bool idle(cuda::GPUPool::instance()->idle());
                    //std::map<void *, int>::iterator j(_device_map.find(map_i->second.address));
                    //running slave thread
                    if (map_i->second.device_id == cuda_get_device() && ! idle)
                    {
                        cuda_download(map_i->second.device, map_i->second.address, map_i->second.bytes);
                    }
                    else
                    {
                        if (! idle)
                        {
                            throw InternalError("MemoryBackend<tags::GPU::CUDA>::download Data is located on another device!");
                        }
                        //running main thread -> switch to slave
                        else
                        {
                            cuda::DownloadTask dt(map_i->second.device, map_i->second.address, map_i->second.bytes);
                            cuda::GPUPool::instance()->enqueue(dt, map_i->second.device_id).wait();
                        }
                    }
                }
            }
        }

        void * alloc(void * memid, void * address, unsigned long bytes)
        {
            std::pair<std::multimap<void *, void *>::iterator, std::multimap<void *, void *>::iterator> range(_address_map.equal_range(address));
            if (range.first == range.second)
            //address not known
            //if (i == _address_map.end())
            {
                void * device(0);
                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (! idle)
                {
                    device = cuda_malloc(bytes);
                }
                //running in main thread -> switch to slave thread
                else
                {
                    cuda::AllocTask at(&device, bytes);
                    cuda::GPUPool::instance()->enqueue(at, cuda_get_device()).wait();
                }
                if (device == 0)
                    throw InternalError("MemoryBackend<tags::GPU::CUDA>::alloc CudaMallocError!");
                _id_map.insert(std::pair<void *, Chunk>(memid, Chunk(address, device, bytes, cuda_get_device())));
                _address_map.insert(std::pair<void *, void *>(address, device));
                _device_map.insert(std::pair<void *, int>(address, cuda_get_device()));
                return device;
            }
            //address already known on some device
            else
            {
                std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
                for (std::multimap<void *, Chunk>::iterator map_i(range.first) ; map_i != range.second ; ++map_i)
                {
                    if (map_i->second.address == address && map_i->second.device_id == cuda_get_device())
                        return map_i->second.device;
                }
                //address was not known on our device, allocate anyway on the actual device
                {
                    void * device(0);
                    bool idle(cuda::GPUPool::instance()->idle());
                    //running in slave thread
                    if (! idle)
                    {
                        device = cuda_malloc(bytes);
                    }
                    //running in main thread -> switch to slave thread
                    else
                    {
                        cuda::AllocTask at(&device, bytes);
                        cuda::GPUPool::instance()->enqueue(at, cuda_get_device()).wait();
                    }
                    if (device == 0)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::alloc CudaMallocError!");
                    _id_map.insert(std::pair<void *, Chunk>(memid, Chunk(address, device, bytes, cuda_get_device())));
                    _address_map.insert(std::pair<void *, void *>(address, device));
                    _device_map.insert(std::pair<void *, int>(address, cuda_get_device()));
                    return device;
                }
            }
        }

        void free(void * memid)
        {
            std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
            bool idle(cuda::GPUPool::instance()->idle());
            for (std::multimap<void *, Chunk>::iterator i(range.first) ; i != range.second ; ++i)
            {
                //running in slave thread
                if (i->second.device_id == cuda_get_device() && ! idle)
                {
                    cuda_free(i->second.device);
                }
                else
                {
                    if (! idle)
                    {
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::free Data is located on another device!");
                    }
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::FreeTask ft(i->second.device);
                        cuda::GPUPool::instance()->enqueue(ft, i->second.device_id).wait();
                    }
                }
                _address_map.erase(i->second.address);
                _device_map.erase(i->second.address);
            }
            _id_map.erase(range.first, range.second);
        }

        void copy(void * /*src_id*/, void * src_address, void * dest_id,
                void * dest_address, unsigned long bytes)
        {
            std::map<void *, void *>::iterator src_i(_address_map.find(src_address));
            std::map<void *, void *>::iterator dest_i(_address_map.find(dest_address));
            if (src_i == _address_map.end() || dest_i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy address not found!");
            }
            else
            {
                std::map<void *, int>::iterator j_source(_device_map.find(src_address));
                std::map<void *, int>::iterator j_dest(_device_map.find(dest_address));
                void * dest_device (0);
                if (j_source->second != j_dest->second)
                {
                    //throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy src and dest on different devices!");
                    this->free(dest_id);
                    cuda::AllocTask at(&dest_device, bytes);
                    cuda::GPUPool::instance()->enqueue(at, j_source->second).wait();
                    if (dest_device == 0)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::alloc CudaMallocError!");
                    _id_map.insert(std::pair<void *, Chunk>(dest_id, Chunk(dest_address, dest_device, bytes, j_source->second)));
                    _address_map.insert(std::pair<void *, void *>(dest_address, dest_device));
                    _device_map.insert(std::pair<void *, int>(dest_address, j_source->second));
                }
                else
                    dest_device = dest_i->second;

                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j_source->second == cuda_get_device() && ! idle)
                {
                    cuda_copy(src_i->second, dest_device, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::CopyTask ct(src_i->second, dest_device, bytes);
                        cuda::GPUPool::instance()->enqueue(ct, j_source->second).wait();
                    }
                }
            }
        }

        void convert_float_double(void * /*src_id*/, void * src_address, void * /*dest_id*/,
                void * dest_address, unsigned long bytes)
        {
            std::map<void *, void *>::iterator src_i(_address_map.find(src_address));
            std::map<void *, void *>::iterator dest_i(_address_map.find(dest_address));
            if (src_i == _address_map.end() || dest_i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy address not found!");
            }
            else
            {
                std::map<void *, int>::iterator j_source(_device_map.find(src_address));
                std::map<void *, int>::iterator j_dest(_device_map.find(dest_address));
                if (j_source->second != j_dest->second)
                    throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy src and dest on different devices!");

                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j_source->second == cuda_get_device() && ! idle)
                {
                    cuda_convert_float_double(src_i->second, dest_i->second, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::ConvertFloatDoubleTask ct(src_i->second, dest_i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ct, j_source->second).wait();
                    }
                }
            }
        }

        void convert_double_float(void * /*src_id*/, void * src_address, void * /*dest_id*/,
                void * dest_address, unsigned long bytes)
        {
            std::map<void *, void *>::iterator src_i(_address_map.find(src_address));
            std::map<void *, void *>::iterator dest_i(_address_map.find(dest_address));
            if (src_i == _address_map.end() || dest_i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy address not found!");
            }
            else
            {
                std::map<void *, int>::iterator j_source(_device_map.find(src_address));
                std::map<void *, int>::iterator j_dest(_device_map.find(dest_address));
                if (j_source->second != j_dest->second)
                    throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy src and dest on different devices!");

                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j_source->second == cuda_get_device() && ! idle)
                {
                    cuda_convert_double_float(src_i->second, dest_i->second, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::copy Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::ConvertDoubleFloatTask ct(src_i->second, dest_i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ct, j_source->second).wait();
                    }
                }
            }
        }

        void fill(void * /*memid*/, void * address, unsigned long bytes, float proto)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            if (i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::fill address not found!");
            }
            else
            {
                if (proto != 0)
                    throw InternalError("CUDA fill != zero not supported yet!");

                std::map<void *, int>::iterator j(_device_map.find(address));
                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j->second == cuda_get_device() && ! idle)
                {
                    cuda_fill_zero(i->second, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::fill Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::FillTask ft(i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ft, j->second).wait();
                    }
                }
            }
        }

        void fill(void * /*memid*/, void * address, unsigned long bytes, double proto)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            if (i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::GPU::CUDA>::fill address not found!");
            }
            else
            {
                if (proto != 0)
                    throw InternalError("CUDA fill != zero not supported yet!");

                std::map<void *, int>::iterator j(_device_map.find(address));
                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j->second == cuda_get_device() && ! idle)
                {
                    cuda_fill_zero(i->second, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::GPU::CUDA>::fill Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::FillTask ft(i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ft, j->second).wait();
                    }
                }
            }
        }

        bool knows(void * /*memid*/, void * address)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            return (i != _address_map.end());
        }
    };

    MemoryBackend<tags::GPU::CUDA>::MemoryBackend() :
        PrivateImplementationPattern<MemoryBackend<tags::GPU::CUDA>, Single>(new Implementation<MemoryBackend<tags::GPU::CUDA> >)
    {
        CONTEXT("When creating MemoryBackend<tags::GPU::CUDA>:");
    }

    MemoryBackend<tags::GPU::CUDA>::~MemoryBackend()
    {
        CONTEXT("When destroying MemoryBackend<tags::GPU::CUDA>:");
    }

    void * MemoryBackend<tags::GPU::CUDA>::upload(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When uploading data (CUDA):");
        return _imp->upload(memid, address, bytes);
    }

    void MemoryBackend<tags::GPU::CUDA>::download(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When downloading data (CUDA):");
        _imp->download(memid, address, bytes);
    }

    void * MemoryBackend<tags::GPU::CUDA>::alloc(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When allocating data (CUDA):");
        return _imp->alloc(memid, address, bytes);
    }

    void MemoryBackend<tags::GPU::CUDA>::free(void * memid)
    {
        CONTEXT("When freeing data (CUDA):");
        _imp->free(memid);
    }

    void MemoryBackend<tags::GPU::CUDA>::copy(void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        CONTEXT("When copying data (CUDA):");
        _imp->copy(src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryBackend<tags::GPU::CUDA>::convert_float_double(void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        CONTEXT("When converting data (CUDA):");
        _imp->convert_float_double(src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryBackend<tags::GPU::CUDA>::convert_double_float(void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        CONTEXT("When converting data (CUDA):");
        _imp->convert_double_float(src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryBackend<tags::GPU::CUDA>::fill(void * memid, void * address, unsigned long bytes, float proto)
    {
        CONTEXT("When filling data (CUDA):");
        _imp->fill(memid, address, bytes, proto);
    }

    void MemoryBackend<tags::GPU::CUDA>::fill(void * memid, void * address, unsigned long bytes, double proto)
    {
        CONTEXT("When filling data (CUDA):");
        _imp->fill(memid, address, bytes, proto);
    }

    bool MemoryBackend<tags::GPU::CUDA>::knows(void * memid, void * address)
    {
        return _imp->knows(memid, address);
    }

    template class InstantiationPolicy<MemoryBackend<tags::GPU::CUDA>, Singleton>;
    //template class MemoryBackend<tags::GPU::CUDA>;
}

using namespace honei;

static MemoryBackendRegistrator gpu_backend_registrator(tags::GPU::CUDA::memory_value, MemoryBackend<tags::GPU::CUDA>::instance());
