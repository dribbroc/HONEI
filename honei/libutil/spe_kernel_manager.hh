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

#ifndef LIBUTIL_GUARD_SPE_KERNEL_MANAGER_HH
#define LIBUTIL_GUARD_SPE_KERNEL_MANAGER_HH 1

#include <honei/libutil/instantiation_policy.hh>
#include <honei/libutil/spe_kernel.hh>
#include <iostream>

namespace honei
{
    /**
     * SPEKernelManager
     *
     * \ingroup grpcell
     */
    class SPEKernelManager :
        public InstantiationPolicy<SPEKernelManager, Singleton>
    {
        private:
            struct Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// \name Basic Operations
            /// \{

            /// Constructor.
            SPEKernelManager();

            /// Destructor.
            ~SPEKernelManager();

            /// \}

        public:
            friend class InstantiationPolicy<SPEKernelManager, Singleton>;

            /// Register a kernel with us.
            void register_kernel(const SPEKernel::Info & info);

            /// \name List of kernels
            /// \{

            typedef libwrapiter::ForwardIterator<SPEKernelManager, SPEKernel::Info> ListIterator;

            /**
             * Returns an iterator pointing to the first kernel in the list.
             */
            ListIterator begin() const;

            /**
             * Returns an iterator pointing behind the last kernel in the list.
             */
            ListIterator end() const;

            /**
             * Returns the number of listed kernels.
             */
            unsigned long size() const;

            /**
             * Return an iterator pointing to the kernel with a given name.
             *
             * \param name Name of the kernel to be found.
             */
            ListIterator find(const std::string & name) const;

            /// \}

            /// \name Map of kernel capabilities
            /// \{

            typedef libwrapiter::ForwardIterator<SPEKernelManager, SPEKernel::Info> MapIterator;

            /**
             * Returns an iterator pointing to the first kernel in a list
             * of kernels that all support a given opcode.
             *
             * \param opcode Opcode that shall be supported.
             */
            MapIterator find(cell::OpCode opcode) const;

            /**
             * Returns an iterator pointing behind the last kernel in a list
             * of kernels that all support a given opcode.
             *
             * \param opcode Opcode that shall be supported.
             */
            MapIterator upper_bound(cell::OpCode opcode) const;

            /// \}
    };

    /**
     * SPEKernelRegistrator
     *
     * \ingroup grpcell
     */
    class SPEKernelRegistrator
    {
        public:
            /**
             * Constructor.
             *
             * \param info The information concerning the to-be-registered kernel.
             *
             * \ingroup grpcell
             */
            SPEKernelRegistrator(const SPEKernel::Info & info)
            {
                // This needs to inline to work around a cyclics dependency of
                // libutil and libcell.

                SPEKernelManager::instance()->register_kernel(info);
            }
    };
}

#endif
