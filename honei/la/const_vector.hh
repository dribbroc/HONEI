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

#include <honei/util/tags.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/stringify.hh>
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

            /**
             * \brief ConstVector::ConstVectorElementIterator is a simple const iterator implementation for const vectors.
             *
             * \ingroup grpvector
             **/
            class ConstVectorElementIterator
            {
                private:
                    /// Our parent vector.
                    const ConstVector<DataType_> & _vector;

                    /// Our index.
                    unsigned long _index;

                public:
                    /// \name Constructors
                    /// \{

                    /**
                     * Constructor.
                     *
                     * \param vector The parent vector that is referenced by the iterator.
                     * \param index The index into the vector.
                     **/
                    ConstVectorElementIterator(const ConstVector<DataType_> & vector, unsigned long index) :
                        _vector(vector),
                        _index(index)
                    {
                    }

                    /// Copy-constructor.
                    ConstVectorElementIterator(const ConstVector<DataType_>::ConstVectorElementIterator & other) :
                        _vector(other._vector),
                        _index(other._index)
                    {
                    }

                    /// Destructor.
                    ~ConstVectorElementIterator()
                    {
                    }

                    /// \}

                    /// \name Forward iterator interface
                    /// \{

                    /// Preincrement operator.
                    ConstVectorElementIterator & operator++ ()
                    {
                        CONTEXT("When incrementing iterator by one:");

                        ++_index;

                        return *this;
                    }

                    /// In-place-add operator.
                    ConstVectorElementIterator & operator+= (const unsigned long step)
                    {
                        CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                        _index += step;

                        return *this;
                    }

                    /// Dereference operator that returns an unassignable reference.
                    const DataType_ & operator* () const
                    {
                        CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                        return _vector[_index];
                    }

                    /// Less-than operator.
                    bool operator< (const ConstVectorElementIterator & other) const
                    {
                        return _index < other.index();
                    }

                    /// Equality operator.
                    bool operator== (const ConstVectorElementIterator & other) const
                    {
                        return ((&_vector == other.parent()) && (_index == other.index()));
                    }

                    /// Inequality operator.
                    bool operator!= (const ConstVectorElementIterator & other) const
                    {
                        return ((&_vector != other.parent()) || (_index != other.index()));
                    }

                    /// \}

                    /// \name IteratorTraits interface
                    /// \{

                    /// Returns our index.
                    virtual unsigned long index() const
                    {
                        return _index;
                    }

                    /// Returns a pointer to our parent container.
                    virtual const ConstVector<DataType_> * parent() const
                    {
                        return &_vector;
                    }
                    /// \}
            };

            /// Returns const iterator pointing to the first element of the vector.
            ConstVectorElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            ConstVectorElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the vector.
            ConstVectorElementIterator element_at(unsigned long index) const;

            /// Returns our size.
            unsigned long size() const;

            /// Retrieves element by index, zero-based, unassignable.
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

    /**
     * Equality operator for ConstVector.
     *
     * Compares if corresponding elements of two const vectors are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const ConstVector<DataType_> & a, const ConstVector<DataType_> & b);

    /**
     * Output operator for ConstVector.
     *
     * Outputs a ConstVector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const ConstVector<DataType_> & vector);

    extern template class ConstVector<float>;

    extern template class ConstVector<double>;
}
#endif
