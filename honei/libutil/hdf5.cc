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

#include <honei/libutil/hdf5.hh>
#include <honei/libutil/instantiation_policy-impl.hh>
#include <honei/libutil/private_implementation_pattern-impl.hh>

#include <map>
#include <tr1/memory>

namespace honei
{
    HDF5Error::HDF5Error(const std::string & function, int reason) :
        ExternalError("libhdf5", function + " failed, returned " + stringify(reason))
    {
    }

    RankExceededError::RankExceededError(unsigned rank) :
        Exception("HDF5SimpleDataSpace created with rank '" + stringify(rank) + "' exceeded its rank")
    {
    }

    template <> struct Implementation<HDF5File> :
        public InstantiationPolicy<Implementation<HDF5File>, NonCopyable>
    {
        typedef std::map<std::string, std::tr1::shared_ptr<Implementation<HDF5Group> > > GroupMap;

        /// Our HDF5 file id.
        const hid_t id;

        /// Our file name.
        const std::string name;

        /// Our access mode.
        const int mode;

        /// Our map of group names to groups.
        GroupMap groups;

        /// Constructor.
        Implementation(const std::string & n, int m) :
            id(H5Fcreate(n.c_str(), m, H5P_DEFAULT, H5P_DEFAULT)),
            name(n),
            mode(m)
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Fcreate", id);
        }

        /// Constructor.
        Implementation(const std::string & n) :
            id(H5Fopen(n.c_str(), hdf5am_read_write, H5P_DEFAULT)),
            name(n),
            mode(hdf5am_read_write)
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Fopen", id);
        }

        /// Destructor.
        ~Implementation()
        {
            herr_t result(H5Fclose(id));

            if (0 > result)
                throw HDF5Error("H5Fclose", result);
        }
    };

    HDF5File::HDF5File(const std::string & name) :
        PrivateImplementationPattern<HDF5File, Shared>(new Implementation<HDF5File>(name))
    {
    }

    HDF5File::HDF5File(Implementation<HDF5File> * imp) :
        PrivateImplementationPattern<HDF5File, Shared>(imp)
    {
    }

    HDF5File
    HDF5File::create(const std::string & name, int mode)
    {
        return HDF5File(new Implementation<HDF5File>(name, mode));
    }


    HDF5Group
    HDF5File::group(const std::string & name)
    {
        Implementation<HDF5File>::GroupMap::iterator g(_imp->groups.find(name));

        if (_imp->groups.end() != g)
        {
            return HDF5Group(g->second.get());
        }
        else
        {
            HDF5Group result(*this, name);

            _imp->groups.insert(std::make_pair(name, result._imp));

            return result;
        }
    }

    hid_t
    HDF5File::id() const
    {
        return _imp->id;
    }

    template <> struct Implementation<HDF5Group> :
        public InstantiationPolicy<Implementation<HDF5Group>, NonCopyable>
    {
        /// Our HDF5 file.
        const HDF5File file;

        /// Our name.
        const std::string name;

        /// Our HDF5 id.
        const hid_t id;

        /**
         * Constructor.
         *
         * Opens a named group in a given HDF5File.
         */
        Implementation(const HDF5File & f, const std::string & n) :
            file(f),
            name(n),
            id(H5Gopen(file.id(), n.c_str()))
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Gopen", id);
        }

        /**
         * Constructor.
         *
         * Creates a named group in a given HDF5File.
         */
        Implementation(const HDF5File & f, const std::string & n, std::size_t unused) :
            file(f),
            name(n),
            id(H5Gcreate(file.id(), n.c_str(), 0))
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Gcreate", id);
        }

        /// Destructor.
        ~Implementation()
        {
            herr_t result(H5Gclose(id));

            if (0 > result)
                throw HDF5Error("H5Gclose", result);
        }
    };

    HDF5Group::HDF5Group(const HDF5File & file, const std::string & name) :
        PrivateImplementationPattern<HDF5Group, Shared>(new Implementation<HDF5Group>(file, name))
    {
    }

    HDF5Group::HDF5Group(Implementation<HDF5Group> * imp) :
        PrivateImplementationPattern<HDF5Group, Shared>(imp)
    {
    }

    HDF5Group
    HDF5Group::create(const HDF5File & file, const std::string & name)
    {
        return HDF5Group(new Implementation<HDF5Group>(file, name, 0));
    }

    template <> struct Implementation<HDF5LoadedDataSpace>
    {
        /// Our HDF5 data space id.
        hid_t id;

        /// Constructor.
        Implementation(hid_t i) :
            id(i)
        {
        }

        /// Destructor.
        ~Implementation()
        {
            if (H5I_INVALID_HID != id)
            {
                herr_t result(H5Sclose(id));

                if (0 > result)
                    throw HDF5Error("H5Sclose", result);
            }
        }

        /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
        Implementation(const Implementation &);

        /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
        Implementation & operator= (const Implementation &);
    };

    HDF5DataSpace::~HDF5DataSpace()
    {
    }

    HDF5LoadedDataSpace::HDF5LoadedDataSpace(hid_t id) :
        PrivateImplementationPattern<HDF5LoadedDataSpace, Shared>(new Implementation<HDF5LoadedDataSpace>(id))
    {
    }

    HDF5LoadedDataSpace::~HDF5LoadedDataSpace()
    {
    }

    hid_t
    HDF5LoadedDataSpace::id() const
    {
        return _imp->id;
    }

    template <> struct Implementation<HDF5SimpleDataSpace> :
        public Implementation<HDF5LoadedDataSpace>
    {
        /// Our rank.
        const unsigned rank;

        /// Our dimensions.
        hsize_t * const dimensions;

        /// Our current index into the dimensions.
        unsigned current;

        /// Constructor
        Implementation(unsigned r) :
            Implementation<HDF5LoadedDataSpace>(H5I_INVALID_HID),
            rank(r),
            dimensions(new hsize_t[rank]),
            current(0)
        {
            std::fill_n(dimensions, rank, 0);
        }

        /// Destructor.
        ~Implementation()
        {
            delete[] dimensions;
        }
    };

    HDF5SimpleDataSpace::HDF5SimpleDataSpace(unsigned rank) :
        PrivateImplementationPattern<HDF5SimpleDataSpace, Shared>(new Implementation<HDF5SimpleDataSpace>(rank))
    {
    }

    HDF5SimpleDataSpace::~HDF5SimpleDataSpace()
    {
    }

    HDF5SimpleDataSpace &
    HDF5SimpleDataSpace::operator[] (hsize_t dimension)
    {
        if (_imp->rank == _imp->current)
            throw RankExceededError(_imp->rank);

        _imp->dimensions[_imp->current] = dimension;
        ++_imp->current;

        if (_imp->rank == _imp->current)
        {
            _imp->id = H5Screate_simple(_imp->rank, _imp->dimensions, 0);

            if (_imp->id == H5I_INVALID_HID)
                throw HDF5Error("H5Screate_simple", _imp->id);
        }

        return *this;
    }

    hid_t
    HDF5SimpleDataSpace::id() const
    {
        return _imp->id;
    };

    template <> struct Implementation<HDF5DataSetBase> :
        public InstantiationPolicy<HDF5DataSetBase, NonCopyable>
    {
        /// Our HDF5 file.
        const HDF5File file;

        /// Our name.
        const std::string name;

        /// Our id.
        const hid_t id;

        /// Our type id.
        const hid_t type_id;

        /// Our data space.
        const HDF5DataSpace * data_space;

        /// Constructor.
        Implementation(const HDF5File & f, const std::string & n, const HDF5SimpleDataSpace & d, hid_t t) :
            file(f),
            name(n),
            id(H5Dcreate(file.id(), name.c_str(), t, d.id(), H5P_DEFAULT)),
            type_id(t),
            data_space(new HDF5SimpleDataSpace(d))
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Dcreate", id);
        }

        /// Constructor.
        Implementation(const HDF5File & f, const std::string & n) :
            file(f),
            name(n),
            id(H5Dopen(file.id(), name.c_str())),
            type_id(H5Dget_type(id)),
            data_space(new HDF5LoadedDataSpace(H5Dget_space(id)))
        {
            if (H5I_INVALID_HID == id)
                throw HDF5Error("H5Dopen", id);
        }

        /// Destructor.
        ~Implementation()
        {
            herr_t result(H5Dclose(id));

            if (0 > result)
                throw HDF5Error("H5Dclose", result);
        }
    };

    HDF5DataSetBase::HDF5DataSetBase(const HDF5File & file, const std::string & name, const HDF5SimpleDataSpace & data_space,
            hid_t type_id) :
        PrivateImplementationPattern<HDF5DataSetBase, Shared>(new Implementation<HDF5DataSetBase>(file, name, data_space, type_id))
    {
    }

    HDF5DataSetBase::HDF5DataSetBase(const HDF5File & file, const std::string & name) :
        PrivateImplementationPattern<HDF5DataSetBase, Shared>(new Implementation<HDF5DataSetBase>(file, name))
    {
    }

    HDF5DataSetBase::~HDF5DataSetBase()
    {
    }

    hid_t
    HDF5DataSetBase::id() const
    {
        return _imp->id;
    }
}
