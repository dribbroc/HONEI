/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/backends/cell/ppe/spe_kernel.hh>
#include <honei/backends/cell/ppe/spe_kernel_manager.hh>
#include <honei/util/assertion.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <tr1/functional>

namespace
{
    bool is_kernel_name(const std::string & name, const ::honei::SPEKernel::Info & info)
    {
        return name == info.name;
    }
}

namespace honei
{
    template class WrappedForwardIterator<SPEKernelManager::ListIteratorTag, SPEKernel::Info>;

    template class WrappedForwardIterator<SPEKernelManager::MapIteratorTag, SPEKernel::Info>;

    using namespace cell;

    template class InstantiationPolicy<SPEKernelManager, Singleton>;

    template <> struct Implementation<SPEKernelManager>
    {
        typedef std::list<cell::KernelInfo> KernelList;

        typedef std::map<cell::OpCode, KernelList> KernelMapByOpCode;

        /// Our list of kernels.
        KernelList kernels;

        /// Our map of kernel lists for specific instructions.
        KernelMapByOpCode map;

        /// Our mutex.
        Mutex * const mutex;

        /// Constructor.
        Implementation() :
            mutex(new Mutex)
        {
        }

        /// Destructor.
        ~Implementation()
        {
            delete mutex;
        }
    };

    SPEKernelManager::SPEKernelManager() :
        PrivateImplementationPattern<SPEKernelManager, Single>(new Implementation<SPEKernelManager>)
    {
    }

    SPEKernelManager::~SPEKernelManager()
    {
    }

    void
    SPEKernelManager::register_kernel(const SPEKernel::Info & info)
    {
        Lock l(*_imp->mutex);

        _imp->kernels.push_back(info);
        for (const OpCode * c(info.capabilities.opcodes), * c_end(c + info.capabilities.opcode_count) ;
                c != c_end ; ++c)
        {
            _imp->map[*c].push_back(info);
        }
    }

    SPEKernelManager::ListIterator
    SPEKernelManager::begin() const
    {
        Lock l(*_imp->mutex);

        return SPEKernelManager::ListIterator(_imp->kernels.begin());
    }

    SPEKernelManager::ListIterator
    SPEKernelManager::end() const
    {
        Lock l(*_imp->mutex);

        return SPEKernelManager::ListIterator(_imp->kernels.end());
    }

    unsigned long
    SPEKernelManager::size() const
    {
        Lock l(*_imp->mutex);

        return _imp->kernels.size();
    }

    SPEKernelManager::ListIterator
    SPEKernelManager::find(const std::string & name) const
    {
        Lock l(*_imp->mutex);

        return SPEKernelManager::ListIterator(std::find_if(_imp->kernels.begin(), _imp->kernels.end(),
                    std::tr1::bind(is_kernel_name, name, std::tr1::placeholders::_1)));
    }

    SPEKernelManager::MapIterator
    SPEKernelManager::find(cell::OpCode opcode) const
    {
        Lock l(*_imp->mutex);

        return SPEKernelManager::MapIterator(_imp->map[opcode].begin());
    }

    SPEKernelManager::MapIterator
    SPEKernelManager::upper_bound(cell::OpCode opcode) const
    {
        Lock l(*_imp->mutex);

        return MapIterator(_imp->map[opcode].end());
    }
}
