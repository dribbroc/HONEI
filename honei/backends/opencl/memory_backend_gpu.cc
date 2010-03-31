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
#include <honei/backends/opencl/opencl_backend.hh>
#include <honei/util/memory_arbiter.hh>

namespace honei
{
    /**
     * \brief MemoryBackend<tags::OpenCL::GPU>::Implementation is the private implementation class for
     * MemoryBackend<tags::OpenCL::GPU>.
     *
     * \ingroup grpvector
     */
    template <> struct Implementation<MemoryBackend<tags::OpenCL::GPU> >
    {
        struct Chunk
        {
            public:
                Chunk(void * a, cl_mem d, unsigned long b) :
                    address(a),
                    device(d),
                    bytes(b)
            {
            }
                void * address;
                cl_mem device;
                unsigned long bytes;
        };

        std::map<void *, void *> _address_map;

        std::multimap<void *, Chunk> _id_map;

        void * upload(void * memid, void * address, unsigned long bytes)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            //address unknown
            if (i == _address_map.end())
            {
                void * device(this->alloc(memid, address, bytes));
                DCQ dcq = OpenCLBackend::instance()->prepare_device(CL_DEVICE_TYPE_GPU);
                clEnqueueWriteBuffer(dcq.command_queue, (cl_mem)device, CL_TRUE, 0, bytes, address, 0, NULL, NULL);
                return device;
            }
            else
            {
                // address already
                return i->second;
            }
        }

        void download(void * memid, void * address, unsigned long bytes)
        {
            std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
            if (range.first == range.second)
            {
                //address not known
                throw InternalError("MemoryBackend<tags::OpenCL::GPU>::download address not found!");
            }
            else
            {
                for (std::multimap<void *, Chunk>::iterator map_i(range.first) ; map_i != range.second ; ++map_i)
                {
                    std::map<void *, void *>::iterator i(_address_map.find(address));
                    DCQ dcq = OpenCLBackend::instance()->prepare_device(CL_DEVICE_TYPE_GPU);
                    clEnqueueReadBuffer(dcq.command_queue, map_i->second.device, CL_TRUE, 0, map_i->second.bytes, map_i->second.address, 0, NULL, NULL);
                    clFinish(dcq.command_queue);
                }
            }
        }

        void * alloc(void * memid, void * address, unsigned long bytes)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            //address not known
            if (i == _address_map.end())
            {
                cl_mem device(0);
                DCQ dcq = OpenCLBackend::instance()->prepare_device(CL_DEVICE_TYPE_GPU);
                device = OpenCLBackend::instance()->create_buffer(bytes, dcq.context, address);
                if (device == 0)
                    throw InternalError("MemoryBackend<tags::OpenCL::GPU>::alloc CudaMallocError!");
                _id_map.insert(std::pair<void *, Chunk>(memid, Chunk(address, device, bytes)));
                _address_map.insert(std::pair<void *, void *>(address, device));
                return device;
            }
            //address already known on some device
            else
            {
                return i->second;
                //throw InternalError("MemoryBackend<tags::OpenCL::GPU>::alloc address already known!");
            }
        }

        void free(void * memid)
        {
            std::multimap<void *, Chunk>::iterator i;
            std::pair<std::multimap<void *, Chunk>::iterator, std::multimap<void *, Chunk>::iterator> range(_id_map.equal_range(memid));
            for (std::multimap<void *, Chunk>::iterator i(range.first) ; i != range.second ; ++i)
            {
                clReleaseMemObject(i->second.device);
                _address_map.erase(i->second.address);
            }
            _id_map.erase(range.first, range.second);
        }

        void copy(void * src_id, void * src_address, void * dest_id,
                void * dest_address, unsigned long bytes)
        {
            /*std::map<void *, void *>::iterator src_i(_address_map.find(src_address));
            std::map<void *, void *>::iterator dest_i(_address_map.find(dest_address));
            if (src_i == _address_map.end() || dest_i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::OpenCL::GPU>::copy address not found!");
            }
            else
            {
                std::map<void *, int>::iterator j_source(_device_map.find(src_address));
                std::map<void *, int>::iterator j_dest(_device_map.find(dest_address));
                if (j_source->second != j_dest->second)
                    throw InternalError("MemoryBackend<tags::OpenCL::GPU>::copy src and dest on different devices!");

                bool idle(cuda::GPUPool::instance()->idle());
                //running in slave thread
                if (j_source->second == cuda_get_device() && ! idle)
                {
                    cuda_copy(src_i->second, dest_i->second, bytes);
                }
                else
                {
                    if (! idle)
                        throw InternalError("MemoryBackend<tags::OpenCL::GPU>::copy Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::CopyTask ct(src_i->second, dest_i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ct, j_source->second)->wait();
                    }
                }
            }*/
        }

        void fill(void * memid, void * address, unsigned long bytes, float proto)
        {
            /*std::map<void *, void *>::iterator i(_address_map.find(address));
            if (i == _address_map.end())
            {
                throw InternalError("MemoryBackend<tags::OpenCL::GPU>::fill address not found!");
            }
            else
            {
                if (proto != 0)
                    throw InternalError("OpenCL::GPU fill != zero not supported yet!");

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
                        throw InternalError("MemoryBackend<tags::OpenCL::GPU>::fill Data is located on another device!");
                    //running in master thread -> switch to slave thread
                    else
                    {
                        cuda::FillTask ft(i->second, bytes);
                        cuda::GPUPool::instance()->enqueue(ft, j->second)->wait();
                    }
                }
            }*/
        }

        bool knows(void * memid, void * address)
        {
            std::map<void *, void *>::iterator i(_address_map.find(address));
            return (i != _address_map.end());
        }
    };

    MemoryBackend<tags::OpenCL::GPU>::MemoryBackend() :
        PrivateImplementationPattern<MemoryBackend<tags::OpenCL::GPU>, Single>(new Implementation<MemoryBackend<tags::OpenCL::GPU> >)
    {
        CONTEXT("When creating MemoryBackend<tags::OpenCL::GPU>:");
    }

    MemoryBackend<tags::OpenCL::GPU>::~MemoryBackend()
    {
        CONTEXT("When destroying MemoryBackend<tags::OpenCL::GPU>:");
    }

    void * MemoryBackend<tags::OpenCL::GPU>::upload(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When uploading data (OpenCL::GPU):");
        return _imp->upload(memid, address, bytes);
    }

    void MemoryBackend<tags::OpenCL::GPU>::download(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When downloading data (OpenCL::GPU):");
        _imp->download(memid, address, bytes);
    }

    void * MemoryBackend<tags::OpenCL::GPU>::alloc(void * memid, void * address, unsigned long bytes)
    {
        CONTEXT("When allocating data (OpenCL::GPU):");
        return _imp->alloc(memid, address, bytes);
    }

    void MemoryBackend<tags::OpenCL::GPU>::free(void * memid)
    {
        CONTEXT("When freeing data (OpenCL::GPU):");
        _imp->free(memid);
    }

    void MemoryBackend<tags::OpenCL::GPU>::copy(void * src_id, void * src_address, void * dest_id,
            void * dest_address, unsigned long bytes)
    {
        CONTEXT("When copying data (OpenCL::GPU):");
        _imp->copy(src_id, src_address, dest_id, dest_address, bytes);
    }

    void MemoryBackend<tags::OpenCL::GPU>::convert_float_double(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes)
    {
    }

    void MemoryBackend<tags::OpenCL::GPU>::convert_double_float(void * src_id, void * src_address, void * dest_id,
                    void * dest_address, unsigned long bytes)
    {
    }

    void MemoryBackend<tags::OpenCL::GPU>::fill(void * memid, void * address, unsigned long bytes, float proto)
    {
        CONTEXT("When filling data (OpenCL::GPU):");
        _imp->fill(memid, address, bytes, proto);
    }

    bool MemoryBackend<tags::OpenCL::GPU>::knows(void * memid, void * address)
    {
        return _imp->knows(memid, address);
    }

    template class InstantiationPolicy<MemoryBackend<tags::OpenCL::GPU>, Singleton>;
    //template class MemoryBackend<tags::OpenCL::GPU>;
}

using namespace honei;

static MemoryBackendRegistrator opencl_gpu_backend_registrator(tags::OpenCL::GPU::memory_value, MemoryBackend<tags::OpenCL::GPU>::instance());
