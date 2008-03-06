/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_INSTANTIATION_POLICY_HH
#define LIBUTIL_GUARD_INSTANTIATION_POLICY_HH 1

namespace honei
{
    /**
     * \{
     * \name Instantiation policy tags
     *
     * Possible tags for InstantiationPolicy.
     */

    /**
     * Tag for an instantiation policy that does not allow copying.
     */
    struct NonCopyable;

    /**
     * Tag for an instantiation policy that does not allow instantiation at all.
     */
    struct NonInstantiable;

    /**
     * Tag for an instantiation policy that does not allow more than one instance of the class.
     */
    struct Singleton;

    /// \}

    /**
     * \{
     * 
     * InstantiationPolicy is a utility base class that can be used to restrict
     * instantiation behaviour of its descendants.
     */

    template <typename T_, typename Method_> class InstantiationPolicy;

    template <typename T_> class InstantiationPolicy<T_, NonCopyable>
    {
        private:
            /// Unwanted copy constructor: Do not implement!
            InstantiationPolicy(const InstantiationPolicy &);

            /// Unwanted copy assignment operator: Do not implement!
            InstantiationPolicy & operator= (const InstantiationPolicy &);

        public:
            /// Default constructor.
            InstantiationPolicy()
            {
            }
    };

    template <typename T_> class InstantiationPolicy<T_, NonInstantiable>
    {
        private:
            /// Unwanted copy constructor: Do not implement!
            InstantiationPolicy(const InstantiationPolicy &);

            /// Unwanted copy assignment operator: Do not implement!
            InstantiationPolicy & operator= (const InstantiationPolicy &);

            /// Unwanted default constructor: Do not implement!
            InstantiationPolicy();
    };

    template <typename T_> class InstantiationPolicy<T_, Singleton>
    {
        private:
            class DeleteOnDestruction;

            friend class DeleteOnDestruction;

            /// Returns a pointer to our instance pointer.
            static T_ * * _instance_ptr();

            /// Deletes the object of T_ that is pointed at by ptr.
            static void _delete(T_ * const ptr);

            /// Unwanted copy constructor: Do not implement!
            InstantiationPolicy(const InstantiationPolicy &);

            /// Unwanted copy assignment operator: Do not implement!
            InstantiationPolicy & operator= (const InstantiationPolicy &);

        protected:
            /// Default constructor.
            InstantiationPolicy()
            {
            }

        public:
            /// Returns the instance.
            static T_ * instance();

    };

    /// \}
}

#endif
