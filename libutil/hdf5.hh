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

#ifndef LIBUTIL_GUARD_HDF5_HH
#define LIBUTIL_GUARD_HDF5_HH 1

#include <libutil/exception.hh>
#include <libutil/log.hh>
#include <libutil/stringify.hh>

#include <string>
#include <tr1/memory>

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
        static inline hid_t type_id() { return H5T_NATIVE_FLOAT; }
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
    class HDF5File
    {
        private:
            /// Our implementation class.
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

        public:
            /**
             * Constructor.
             *
             * \param name The new file's name.
             * \param mode The new file's creation mode.
             */
            HDF5File(const std::string & name, int mode = 0);

            /// Return our HDF5 object id.
            hid_t id() const;
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
        protected:
            /// Constructor.
            HDF5DataSpace();

        public:
            /// Destructor.
            virtual ~HDF5DataSpace();

            /// Return a copy of this data space.
            virtual HDF5DataSpace * copy() const;

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
        public HDF5DataSpace
    {
        private:
            /// Our implementation class.
            struct Implementation;

            /// Our implementaion.
            std::tr1::shared_ptr<Implementation> _imp;

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

            /// Return a copy of this data space.
            virtual HDF5SimpleDataSpace * copy() const;

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
    class HDF5DataSetBase
    {
        private:
            /// Our implementation class.
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

        public:
            /**
             * Constructor.
             *
             * Create a new HDF5 dataset in a given HDF5 file.
             *
             * \param file File in which the set shall be created.
             * \param name Name of the data set.
             * \param data_space The data space which this set shall use.
             * \param type_id The type id of the set's elements.
             */
            HDF5DataSetBase(const HDF5File & file, const std::string & name, const HDF5DataSpace & data_space,
                    hid_t type_id);

            /// Destructor.
            virtual ~HDF5DataSetBase();

            /// Return our HDF5 object id.
            hid_t id() const;
    };

    /**
     * HDF5DataSet is the wrapper class template for libhdf5's (H5D-prefixed) data set functions.
     *
     * Creates a new HDF5 data set of given data type and data spacee and closes it once no HDF5
     * wrapper references it anymore.
     *
     * \ingroup grphdf5
     */
    template <typename DT_> class HDF5DataSet :
        public HDF5DataSetBase
    {
        public:
            HDF5DataSet(HDF5File & file, const std::string & name, const HDF5DataSpace & data_space) :
                HDF5DataSetBase(file, name, data_space, HDF5Type<DT_>::type_id())
            {
            }

            HDF5DataSet<DT_> & operator<< (DT_ * values)
            {
                herr_t result;

                if (0 > (result = H5Dwrite(id(), HDF5Type<DT_>::type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values)))
                    throw HDF5Error("H5Dwrite", result);

                return *this;
            }
    };
}

#endif
