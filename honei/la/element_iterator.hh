/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
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

#ifndef LIBLA_GUARD_ELEMENT_ITERATOR_HH
#define LIBLA_GUARD_ELEMENT_ITERATOR_HH 1

#include <honei/la/banded_matrix-fwd.hh>
#include <honei/la/const_vector-fwd.hh>
#include <honei/la/dense_matrix-fwd.hh>
#include <honei/la/dense_matrix_tile-fwd.hh>
#include <honei/la/dense_vector-fwd.hh>
#include <honei/la/dense_vector_range-fwd.hh>
#include <honei/la/dense_vector_slice-fwd.hh>
#include <honei/la/sparse_vector-fwd.hh>
#include <honei/la/sparse_matrix-fwd.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/shared_array-fwd.hh>

#include <iterator>
#include <tr1/memory>

namespace honei
{
    // New implementation
    namespace container
    {
        struct Matrix;

        struct MatrixTile;

        struct Vector;
    }

    namespace storage
    {
        struct Banded;

        struct Const;

        struct Dense;

        struct Sparse;

        struct SparseNonZero;
    }

    template <typename Storage_, typename Container_, typename DataType_> class ElementIterator;

    template <typename Storage_, typename Container_, typename DataType_> class ConstElementIterator;

    template <typename DataType_> class ElementIterator<storage::Banded, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Banded, container::Matrix, DataType_>, Single>
    {
        private:
           /// Constructor.
            ElementIterator(const BandedMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ElementIterator
            /// \{

            friend class ConstElementIterator<storage::Banded, container::Matrix, DataType_>;
            friend class BandedMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ElementIterator<storage::Banded, container::Matrix, float>;

    extern template class ElementIterator<storage::Banded, container::Matrix, double>;

    template <typename DataType_> class ElementIterator<storage::Sparse, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(SparseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ElementIterator
            /// \{

            friend class ConstElementIterator<storage::Sparse, container::Matrix, DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ElementIterator<storage::Sparse, container::Matrix, float>;

    extern template class ElementIterator<storage::Sparse, container::Matrix, double>;

    template <typename DataType_> class ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(SparseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ElementIterator
            /// \{

            friend class ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ElementIterator<storage::SparseNonZero, container::Matrix, float>;

    extern template class ElementIterator<storage::SparseNonZero, container::Matrix, double>;

    template <typename DataType_> class ElementIterator<storage::Dense, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Dense, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(const DenseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ElementIterator
            /// \{

            friend class ConstElementIterator<storage::Dense, container::Matrix, DataType_>;
            friend class DenseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ElementIterator<storage::Dense, container::Matrix, float>;

    extern template class ElementIterator<storage::Dense, container::Matrix, double>;

    template <typename DataType_> class ElementIterator<storage::Dense, container::MatrixTile, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(DenseMatrixTile<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ElementIterator
            /// \{

            friend class ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>;
            friend class DenseMatrixTile<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name MatrixTile element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ElementIterator<storage::Dense, container::MatrixTile, float>;

    extern template class ElementIterator<storage::Dense, container::MatrixTile, double>;

    template <typename DataType_> class ElementIterator<storage::Dense, container::Vector, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Dense, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(const SharedArray<DataType_> & elements, unsigned long index, unsigned long offset,
                    unsigned long stepsize);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class ConstElementIterator<storage::Dense, container::Vector, DataType_>;
            friend class DenseVector<DataType_>;
            friend class DenseVectorRange<DataType_>;
            friend class DenseVectorSlice<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ElementIterator<storage::Dense, container::Vector, float>;

    extern template class ElementIterator<storage::Dense, container::Vector, double>;

    template <typename DataType_> class ElementIterator<storage::Sparse, container::Vector, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(SparseVector<DataType_> & vector, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class ConstElementIterator<storage::Sparse, container::Vector, DataType_>;
            friend class SparseVector<DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Dereference operator returning an unassignable reference.
            //const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ElementIterator<storage::Sparse, container::Vector, float>;

    extern template class ElementIterator<storage::Sparse, container::Vector, double>;

    template <typename DataType_> class ElementIterator<storage::SparseNonZero, container::Vector, DataType_> :
        public PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ElementIterator(SparseVector<DataType_> & vector, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>;
            friend class SparseVector<DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ElementIterator(const ElementIterator &);

            /// Destructor.
            ~ElementIterator();

            /// Assignment operator.
            ElementIterator & operator= (const ElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ElementIterator & operator++ ();

            /// In-place-add operator.
            ElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an assignable reference.
            DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ElementIterator<storage::SparseNonZero, container::Vector, float>;

    extern template class ElementIterator<storage::SparseNonZero, container::Vector, double>;

    template <typename DataType_> class ConstElementIterator<storage::Banded, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Banded, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const BandedMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class BandedMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Banded, container::Matrix, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Banded, container::Matrix, float>;

    extern template class ConstElementIterator<storage::Banded, container::Matrix, double>;

    template <typename DataType_> class ConstElementIterator<storage::Const, container::Vector, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Const, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const ConstVector<DataType_> & vector, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class ConstVector<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Const, container::Vector, float>;

    extern template class ConstElementIterator<storage::Const, container::Vector, double>;

    template <typename DataType_> class ConstElementIterator<storage::Sparse, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Sparse, container::Matrix, float>;

    extern template class ConstElementIterator<storage::Sparse, container::Matrix, double>;

    template <typename DataType_> class ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::SparseNonZero, container::Matrix, float>;

    extern template class ConstElementIterator<storage::SparseNonZero, container::Matrix, double>;

    template <typename DataType_> class ConstElementIterator<storage::Dense, container::Matrix, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Matrix, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const DenseMatrix<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class DenseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Dense, container::Matrix, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Matrix element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Dense, container::Matrix, float>;

    extern template class ConstElementIterator<storage::Dense, container::Matrix, double>;

    template <typename DataType_> class ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const DenseMatrixTile<DataType_> & matrix, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class DenseMatrixTile<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name MatrixTile element iteration interface
            /// \{

            unsigned long column() const;

            unsigned long index() const;

            unsigned long row() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Dense, container::MatrixTile, float>;

    extern template class ConstElementIterator<storage::Dense, container::MatrixTile, double>;

    template <typename DataType_> class ConstElementIterator<storage::Dense, container::Vector, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const SharedArray<DataType_> & elements, unsigned long index, unsigned long offset,
                    unsigned long stepsize);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class DenseVector<DataType_>;
            friend class DenseVectorRange<DataType_>;
            friend class DenseVectorSlice<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Dense, container::Vector, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Dense, container::Vector, float>;

    extern template class ConstElementIterator<storage::Dense, container::Vector, double>;

    template <typename DataType_> class ConstElementIterator<storage::Sparse, container::Vector, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const SparseVector<DataType_> & vector, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class SparseVector<DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::Sparse, container::Vector, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };

    extern template class ConstElementIterator<storage::Sparse, container::Vector, float>;

    extern template class ConstElementIterator<storage::Sparse, container::Vector, double>;

    template <typename DataType_> class ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> :
        public PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstElementIterator(const SparseVector<DataType_> & vector, unsigned long index);

        public:
            /// \name Friends of ConstElementIterator
            /// \{

            friend class SparseVector<DataType_>;
            friend class SparseMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstElementIterator(const ConstElementIterator & other);

            /// Constructor from ElementIterator.
            ConstElementIterator(const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other);

            /// Destructor.
            ~ConstElementIterator();

            /// Assignment operator.
            ConstElementIterator & operator= (const ConstElementIterator &);

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstElementIterator & operator++ ();

            /// In-place-add operator.
            ConstElementIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            const DataType_ & operator* () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstElementIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstElementIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstElementIterator & other) const;

            /// \}

            /// \name Vector element iteration interface
            /// \{

            unsigned long index() const;

            /// \}
    };
}

#endif
