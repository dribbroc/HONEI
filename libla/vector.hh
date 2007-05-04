/* vim: set sw=4 sts=4 et nofoldenable : */

#ifndef BASISKOMPONENTEN_GUARD_TYPES_HH
#define BASISKOMPONENTEN_GUARD_TYPES_HH 1

#include <iterator>
#include <libwrapiter/libwrapiter_forward_iterator.hh>
#include <tr1/memory>

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A Vector is the abstract baseclass for all vector-like types used.
     **/
    template <typename DataType_> class Vector
    {
        public:
            /// Returns iterator pointing to the first element of the vector.
            virtual libwrapiter::ForwardIterator<Vector<DataType_>, DataType_> begin_elements() const = 0;

            /// Returns iterator pointing behind the last element of the vector.
            virtual libwrapiter::ForwardIterator<Vector<DataType_>, DataType_> end_elements() const = 0;

            /// Returns our size.
            virtual unsigned long size() const = 0;

            /// Retrieves element by index, zero-based.
            virtual DataType_ & operator[] (unsigned long index) = 0;
    };

    /**
     * A DenseVector is a vector with O(size**2) elements which keeps its data
     * sequential.
     **/
    template <typename DataType_> class DenseVector :
        public Vector<DataType_>
    {
        private:
            /// Pointer to our elements.
            DataType_ *_elements;

            /// Our size.
            unsigned long _size;

        public:
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Constructor.
            DenseVector(unsigned long size) :
                _elements(new DataType_[size]),
                _size(size)
            {
            }

            /// Destructor.
            ~DenseVector()
            {
                /// \todo shared array. Do not use shared_ptr due to delete/delete[] mismatch.
                delete[] _elements;
            }

            /// Returns iterator pointing to the first element of the vector.
            virtual libwrapiter::ForwardIterator<Vector<DataType_>, DataType_> begin_elements() const
            {
                return libwrapiter::ForwardIterator<Vector<DataType_>, DataType_>(ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual libwrapiter::ForwardIterator<Vector<DataType_>, DataType_> end_elements() const
            {
                return libwrapiter::ForwardIterator<Vector<DataType_>, DataType_>(ElementIteratorImpl<DataType_>(*this, this->size()));
            }

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _size;
            }

            /// Retrieves element by index, zero-based.
            virtual DataType_ & operator[] (unsigned long index)
            {
                return _elements[index];
            }
    };

    /**
     * A DenseVector::ElementIteratorImpl is a simple iterator implementation for dense vectors.
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::ElementIteratorImpl<DataType_> :
        std::iterator<std::forward_iterator_tag, DataType_>
    {
        private:
            const DenseVector<DataType_> & _vector;
            unsigned long _pos;

        public:
            /// Constructor.
            ElementIteratorImpl(const DenseVector<DataType_> & vector, unsigned long pos) :
                _vector(vector),
                _pos(pos)
            {
            }

            /// Constructor.
            ElementIteratorImpl(ElementIteratorImpl<DataType_> const & other) :
                _vector(other._vector),
                _pos(other._pos)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_pos;
                return *this; /// \todo friend for DenseVector!
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<DataType_> operator++ (int)
            {
                ElementIteratorImpl<DataType_> result(*this);
                ++_pos;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImpl<DataType_> & other) const
            {
                return (&_vector == &other._vector) && (_pos == other._pos);
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImpl<DataType_> & other) const
            {
                return ((&_vector != &other._vector) || (_pos != other._pos));
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return _vector._elements[_pos];
            }
    };
}
