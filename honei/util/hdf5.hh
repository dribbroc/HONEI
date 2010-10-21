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

#pragma once
#ifndef LIBUTIL_GUARD_HDF5_HH
#define LIBUTIL_GUARD_HDF5_HH 1

#include <honei/util/exception.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>

#include <string>

#include <hdf5.h>

/**
 * \def H5OPEN
 *
 * Get rid of the abomination that is H5OPEN in what otherwise would be constant expressions.
 *
 * \ingroup grphdf5
 */
#undef H5OPEN
#define H5OPEN

/**
 * \def H5CHECK
 *
 * Get rid of the abomination that is H5CHECK in what otherwise would be constant expressions.
 *
 * \ingroup grphdf5
 */
#undef H5CHECK
#define H5CHECK

namespace honei
{
    // Forward declaration
    class HDF5Group;

    /**
     * HDF5Error is thrown by libhdf5 wrapper classes whenever an error occurs in interfacing libhdf5.
     *
     * \ingroup grpexceptions
     * \ingroup grphdf5
     */
    class HDF5Error :
        public ExternalError
    {
        public:
            /**
             * Constructor.
             *
             * \param function Name of the function that failed.
             * \param reason The reason for the error.
             */
            HDF5Error(const std::string & function, int reason);
    };

    /**
     * RankExceededError is thrown by HDF5SimpleDataSpace when it is instructed to create
     * more dimensions than initially told.
     *
     * \ingroup grpexceptions
     * \ingroup grphdf5
     */
    class RankExceededError :
        public Exception
    {
        public:
            /**
             * Constructor.
             *
             * \param rank Rank with which HDF5SimpleDataSpace was constructed.
             */
            RankExceededError(unsigned rank);
    };

    /**
     * \{
     *
     * HDF5Type is a utility class that allows easy conversion of C++ types to libhd5 type ids.
     *
     * \ingroup grphdf5
     */
    template <typename DT_> struct HDF5Type
    {
    };

    template <> struct HDF5Type<float>
    {
        static inline hid_t storage_type_id() { return H5T_IEEE_F32BE; }
        static inline hid_t memory_type_id() { return H5T_NATIVE_FLOAT; }
    };

    /// \}

    /**
     * Access modes for HDF5 files.
     *
     * \ingroup grphdf5
     */
    enum HDF5AccessMode
    {
        hdf5am_read_only = H5F_ACC_RDONLY,
        hdf5am_read_write = H5F_ACC_RDWR,
        hdf5am_truncate = H5F_ACC_TRUNC,
        hdf5am_exclusive = H5F_ACC_EXCL,
        hdf5am_create = H5F_ACC_CREAT
    };

    /**
     * HDF5File is the wrapper class for libhdf5's (H5F-prefixed) file access functions.
     *
     * Creates a new HDF5 file upon construction and closes it once no HDF5 wrapper
     * references the file anymore.
     *
     * \ingroup grphdf5
     */
    class HDF5File :
        public PrivateImplementationPattern<HDF5File, Shared>
    {
        private:
            /**
             * Constructor.
             *
             * For internal use only.
             */
            HDF5File(Implementation<HDF5File> * imp);

        public:
            friend class HDF5Group;

            ~HDF5File(){}
            /**
             * Constructor.
             *
             * Open an existing HDF5 file of given name.
             *
             * \param name The file's name.
             */
            HDF5File(const std::string & name);

            /**
             * Named Constructor.
             *
             * Create a new HDF5 file of given name.
             *
             * \param name The new file's name.
             * \param mode The new file's creation mode.
             */
            static HDF5File create(const std::string & name, int mode = 0);

            /**
             * Return a named group object.
             *
             * \param name Name of the group whose wrapper shall be returned.
             */
            HDF5Group group(const std::string & name);

            /// Return our HDF5 object id.
            hid_t id() const;
    };

    /**
     * HDF5Group is the wrapper class for libhdf5's (HFG-prefixed) group access functions.
     *
     * Creates or opens a new HDF5 group upon construction and closes it once no HDF5
     * wrapper references the group anymore.
     *
     * \ingroup grphdf5
     */
    class HDF5Group :
        public PrivateImplementationPattern<HDF5Group, Shared>
    {
        private:
            /**
             * Constructor.
             *
             * Opens a given group.
             *
             * \param file The new group's file.
             * \param name The new group's name.
             */
            HDF5Group(const HDF5File & file, const std::string & name);

            /**
             * Constructor.
             *
             * Creates a HDF5Group object from an existing Implementation object.
             *
             * \param imp The group's implementation object.
             */
            HDF5Group(Implementation<HDF5Group> * _imp);

        public:
            friend class HDF5File;

            /// \name Basic operations
            /// \{

            /**
             * Named constructor.
             *
             * Creates a new HDF5 group of given name.
             *
             * \param file The HDF5 file in which the group shall be created.
             * \param name The name of the group that shall be created.
             */
            static HDF5Group create(const HDF5File & file, const std::string & name);

            /// \}
    };

    /**
     * HDF5DataSpace is the base class for wrappers of both simple and complex data spaces
     * as supported by libhdf5.
     *
     * \see HDF5SimpleDataSpace
     *
     * \ingroup grphdf5
     */
    class HDF5DataSpace
    {
        public:
            /// Destructor.
            virtual ~HDF5DataSpace();

            /// Return our HDF5 object id.
            virtual hid_t id() const = 0;
    };

    /**
     * HDF5LoadedDataSpace is the wrapper for existing and loaded data spaces.
     *
     * For use by HDF5 wrapper classes only.
     *
     * \ingroup grphdf5
     */
    class HDF5LoadedDataSpace :
        public HDF5DataSpace,
        public PrivateImplementationPattern<HDF5LoadedDataSpace, Shared>
    {
        public:
            friend class HDF5DataSetBase;

            /**
             * Constructor.
             *
             * Create a HDF5 data space wrapper object from an HDF5 data space id.
             *
             * \param id The data space's object id.
             */
            HDF5LoadedDataSpace(hid_t id);

            /// Destructor.
            virtual ~HDF5LoadedDataSpace();

            /// Return our HDF5 object id.
            virtual hid_t id() const;
    };

    /**
     * HDF5SimpleDataSpace is the wrapper class for libhdf5's (H5S-prefixed) data space functions.
     *
     * Creates a new simple HDF5 data space and closes it once no HDF5 wrapper references
     * it anymore.
     *
     * Due to constraints imposed by the C++ language standard (C++-03), the contruction of a rank \f$r\f$
     * dataspace of dimensions \f$d_0, \dots, d_{r - 1}\f$ works as follows:
     * \code
     * // Let r be 2, d0 be 10 and d1 be 5:
     * HDF5SimpleDataSpace dataspace(r);
     * dataspace[10][5];
     * \endcode
     *
     * \ingroup grphdf5
     */
    class HDF5SimpleDataSpace :
        public HDF5DataSpace,
        public PrivateImplementationPattern<HDF5SimpleDataSpace, Shared>
    {
        public:
            /**
             * Constructor.
             *
             * Create a simple HDF5 data space of given rank.
             *
             * \param rank Rank of the data space.
             */
            HDF5SimpleDataSpace(unsigned rank);

            /// Destructor.
            virtual ~HDF5SimpleDataSpace();

            /**
             * Subscript operator.
             *
             * Add another dimension to the dataspace.
             *
             * \param dimension Size of a dimension that is to be added to the data space.
             */
            HDF5SimpleDataSpace & operator[] (hsize_t dimension);

            /// Return our HDF5 object id.
            virtual hid_t id() const;
    };

    /**
     * HDF5DataSetBase is the base class all data set wrappers.
     *
     * \see HDF5DataSet
     *
     * \ingroup grphdf5
     */
    class HDF5DataSetBase :
        PrivateImplementationPattern<HDF5DataSetBase, Shared>
    {
        protected:
            /**
             * Constructor.
             *
             * Create a new HDF5 data set in a given HDF5 file.
             *
             * \param file The file in which the set shall be created.
             * \param name The name of the data set.
             * \param data_space The data space which this set shall use.
             * \param type_id The type id of the set's elements.
             */
            HDF5DataSetBase(const HDF5File & file, const std::string & name, const HDF5SimpleDataSpace & data_space,
                    hid_t type_id);

            /**
             * Constructor.
             *
             * Open an existing HDF5 data set in a given HDF5 file.
             *
             * \param file The File in which the data set is located.
             * \param name The name of the data set.
             */
            HDF5DataSetBase(const HDF5File & file, const std::string & name);

        public:

            /// Destructor.
            virtual ~HDF5DataSetBase();

            /// Return our HDF5 object id.
            hid_t id() const;
    };

    /**
     * HDF5DataSet is the wrapper class template for libhdf5's (H5D-prefixed) data set functions.
     *
     * Creates a new HDF5 data set of given data type and data space and closes it once no HDF5
     * wrapper references it anymore.
     *
     * \ingroup grphdf5
     */
    template <typename DT_> class HDF5DataSet :
        public HDF5DataSetBase
    {
        private:
            /**
             * Constructor.
             *
             * Create a new data set of given name in a given HDF5 file.
             *
             * \param file The new data set's parent file.
             * \param name The new data set's name.
             * \param data_space The new data set's data space.
             */
            HDF5DataSet(HDF5File & file, const std::string & name, const HDF5SimpleDataSpace & data_space) :
                HDF5DataSetBase(file, name, data_space, HDF5Type<DT_>::storage_type_id())
            {
            }

        public:
            /**
             * Constructor.
             *
             * Open an existing data set of given name from a given HDF5 file.
             *
             * \param file The data set's parent file.
             * \param name The data set's name.
             */
            HDF5DataSet(HDF5File & file, const std::string & name) :
                HDF5DataSetBase(file, name)
            {
            }

            /**
             * Named constructor.
             *
             * Create a new data set of given name in a given HDF5 file.
             *
             * \param file The new data set's parent file.
             * \param name The new data set's name.
             * \param data_space The new data set's data space.
             */
            static inline HDF5DataSet<DT_> create(HDF5File & file, const std::string & name,
                    const HDF5SimpleDataSpace & data_space)
            {
                return HDF5DataSet<DT_>(file, name, data_space);
            }

            /**
             * Serialize data to the data set.
             *
             * \param values The data that shall be written to the data set.
             */
            HDF5DataSet<DT_> & operator<< (const DT_ * values)
            {
                herr_t result;

                if (0 > (result = H5Dwrite(id(), HDF5Type<DT_>::memory_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values)))
                    throw HDF5Error("H5Dwrite", result);

                return *this;
            }

            /**
             * Deserialize data from the data set.
             *
             * \param values The data that shall be read from the data set.
             */
            const HDF5DataSet<DT_> & operator>> (DT_ * values) const
            {
                herr_t result;

                if (0 > (result = H5Dread(id(), HDF5Type<DT_>::memory_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values)))
                    throw HDF5Error("H5Dread", result);

                return *this;
            }
    };
}

#endif
