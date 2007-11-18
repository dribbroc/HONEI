/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libutil/hdf5.hh>

using namespace honei;

HDF5Error::HDF5Error(const std::string & function, int reason) :
    ExternalError("libhdf5", function + " failed, returned " + stringify(reason))
{
}

RankExceededError::RankExceededError(unsigned rank) :
    Exception("HDF5SimpleDataSpace created with rank '" + stringify(rank) + "' exceeded its rank")
{
}

struct HDF5File::Implementation
{
    /// Our libhdf5 file id.
    const hid_t id;

    /// Our file name.
    const std::string name;

    /// Our access mode.
    const int mode;

    /// Constructor.
    Implementation(const std::string & n, int m) :
        id(H5Fcreate(n.c_str(), m, H5P_DEFAULT, H5P_DEFAULT)),
        name(n),
        mode(m)
    {
        if (H5I_INVALID_HID == id)
            throw HDF5Error("H5Fcreate", id);
    }

    /// Destructor.
    ~Implementation()
    {
        herr_t result(H5Fclose(id));

        if (0 > result)
            throw HDF5Error("H5Fclose", result);
    }
};

HDF5File::HDF5File(const std::string & name, int mode) :
    _imp(new Implementation(name, mode))
{
}

hid_t
HDF5File::id() const
{
    return _imp->id;
}

HDF5DataSpace::HDF5DataSpace()
{
}

HDF5DataSpace::~HDF5DataSpace()
{
}

HDF5DataSpace *
HDF5DataSpace::copy() const
{
    return new HDF5DataSpace(*this);
}

hid_t
HDF5DataSpace::id() const
{
    return H5I_INVALID_HID;
}

struct HDF5SimpleDataSpace::Implementation
{
    /// Our rank.
    const unsigned rank;

    /// Our dimensions.
    hsize_t * const dimensions;

    /// Our current index into the dimensions.
    unsigned current;

    /// Our id.
    hid_t id;

    /// Constructor
    Implementation(unsigned r) :
        rank(r),
        dimensions(new hsize_t[rank]),
        current(0)
    {
        std::fill_n(dimensions, rank, 0);
    }

    /// Destructor.
    ~Implementation()
    {
        herr_t result(H5Sclose(id));

        if (0 > result)
            throw HDF5Error("H5Sclose", result);

        delete[] dimensions;
    }
};

HDF5SimpleDataSpace::HDF5SimpleDataSpace(unsigned rank) :
    _imp(new Implementation(rank))
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

HDF5SimpleDataSpace *
HDF5SimpleDataSpace::copy() const
{
    return new HDF5SimpleDataSpace(*this);
}

hid_t
HDF5SimpleDataSpace::id() const
{
    return _imp->id;
}

struct HDF5DataSetBase::Implementation
{
    /// Our HDF5 file.
    const HDF5File file;

    /// Our name.
    const std::string name;

    /// Our data space.
    const HDF5DataSpace * data_space;

    /// Our type id.
    const hid_t type_id;

    /// Our id.
    const hid_t id;

    /// Constructor.
    Implementation(const HDF5File & f, const std::string & n, const HDF5DataSpace & d, hid_t t) :
        file(f),
        name(n),
        data_space(d.copy()),
        type_id(t),
        id(H5Dcreate(file.id(), name.c_str(), type_id, data_space->id(), H5P_DEFAULT))
    {
        if (H5I_INVALID_HID == id)
            throw HDF5Error("H5Dcreate", id);
    }

    /// Destructor.
    ~Implementation()
    {
        delete data_space;
    }
};

HDF5DataSetBase::HDF5DataSetBase(const HDF5File & file, const std::string & name, const HDF5DataSpace & data_space,
        hid_t type_id) :
    _imp(new Implementation(file, name, data_space, type_id))
{
}

HDF5DataSetBase::~HDF5DataSetBase()
{
    herr_t result(H5Dclose(_imp->id));

    if (0 > result)
        throw HDF5Error("H5Dclose", result);
}

hid_t
HDF5DataSetBase::id() const
{
    return _imp->id;
}
