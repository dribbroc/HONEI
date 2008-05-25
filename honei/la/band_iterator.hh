/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef HONEI_GUARD_LA_BAND_ITERATOR_HH
#define HONEI_GUARD_LA_BAND_ITERATOR_HH 1

#include <honei/la/dense_vector-fwd.hh>
#include <honei/la/banded_matrix-fwd.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    namespace type
    {
        struct Banded;
    }

    template <typename Type_, typename DataType_> class BandIterator;

    template <typename Type_, typename DataType_> class ConstBandIterator;

    template <typename DataType_> class BandIterator<type::Banded, DataType_> :
        public PrivateImplementationPattern<BandIterator<type::Banded, DataType_>, Single>
    {
        private:
            /// Constructor.
            BandIterator(Implementation<BandIterator<type::Banded, DataType_> > * imp);

        public:
            /// \name Friends of ConstBandIterator
            /// \{

            friend class ConstBandIterator<type::Banded, DataType_>;
            friend class BandedMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            BandIterator(const BandIterator &);

            /// Destructor.
            ~BandIterator();

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            BandIterator & operator++ ();

            /// In-place-add operator.
            BandIterator & operator+= (const unsigned long step);

            /**
             * Dereference operator returning a DenseVector.
             *
             * \warning Will create the band if it doesn't yet exist.
             */
            DenseVector<DataType_> operator* () const;

            /**
             * Dereference operator that allows access to DenseVector members.
             *
             * \warning Will create the band if it doesn't yet exist.
             */
            DenseVector<DataType_> operator-> () const;

            /// Comparison operator for less-than relations.
            bool operator< (const BandIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const BandIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const BandIterator & other) const;

            /// \}

            /// \name Matrix band iteration interface
            /// \{

            bool exists() const;

            unsigned long index() const;

            /// \}
    };

    extern template class BandIterator<type::Banded, float>;

    extern template class BandIterator<type::Banded, double>;

    template <typename DataType_> class ConstBandIterator<type::Banded, DataType_> :
        public PrivateImplementationPattern<ConstBandIterator<type::Banded, DataType_>, Single>
    {
        private:
            /// Constructor.
            ConstBandIterator(Implementation<ConstBandIterator<type::Banded, DataType_> > * imp);

        public:
            /// \name Friends of ConstBandIterator
            /// \{

            friend class BandedMatrix<DataType_>;

            /// \}

            /// Copy-constructor.
            ConstBandIterator(const ConstBandIterator & other);

            /// Constructor from BandIterator.
            ConstBandIterator(const BandIterator<type::Banded, DataType_> & other);

            /// Destructor.
            ~ConstBandIterator();

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            ConstBandIterator & operator++ ();

            /// In-place-add operator.
            ConstBandIterator & operator+= (const unsigned long step);

            /// Dereference operator returning an unassignable reference.
            /// \todo Change to ConstVector!
            DenseVector<DataType_> operator* () const;

            /// Arrow operator.
            /// \todo Change to ConstVector!
            DenseVector<DataType_> operator-> () const;

            /// Comparison operator for less-than relations.
            bool operator< (const ConstBandIterator & other) const;

            /// Comparison operator for is-equal relations.
            bool operator== (const ConstBandIterator & other) const;

            /// Comparison operator for is-unequal relations.
            bool operator!= (const ConstBandIterator & other) const;

            /// \}

            /// \name Matrix band iteration interface
            /// \{

            bool exists() const;

            unsigned long index() const;

            /// \}
    };

    extern template class ConstBandIterator<type::Banded, float>;

    extern template class ConstBandIterator<type::Banded, double>;
}

#endif
