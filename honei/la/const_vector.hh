/* vim: set sw=4 sts=4 et nofoldenable : */

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

#ifndef LIBLA_GUARD_CONST_VECTOR_HH
#define LIBLA_GUARD_CONST_VECTOR_HH 1

#include <honei/la/vector.hh>
#include <honei/util/tags.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/la/dense_vector.hh>

namespace honei
{
    /**
     * ConstVector is a vector with read only element access.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class ConstVector :
        public PrivateImplementationPattern<ConstVector<DataType_>, Shared>
    {
        private:

        public:
            /// Constructors
            /// \{


            /**
             * Constructor.
             *
             * \param src The source vector, whose data should be used
             */
            ConstVector(const DenseVector<DataType_> & src);


            /// Copy-constructor.
            ConstVector(const ConstVector<DataType_> & other);

            /// Destructor.
            ~ConstVector();

            /// \}


            /*/// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the vector.
            virtual ConstElementIterator element_at(unsigned long index) const;
*/
            /// Returns our size.
            unsigned long size() const;

            /// Retrieves element by index, zero-based, assignable.
            const DataType_ & operator[] (unsigned long index) const;

            /// Return our offset.
            unsigned long offset() const;

            /// Return our memory id
            unsigned long memid() const;

            /// Return the address of our data
            void * address() const;

            /// Request a read lock for our data.
            void * read(tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a read lock for our data.
            void release_read() const;

            /// \}
    };

    extern template class ConstVector<float>;

    extern template class ConstVector<double>;
}

#endif
