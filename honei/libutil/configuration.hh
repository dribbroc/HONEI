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

#ifndef LIBUTIL_GUARD_CONFIGURATION_HH
#define LIBUTIL_GUARD_CONFIGURATION_HH 1

#include <honei/libutil/exception.hh>

#include <string>

namespace honei
{
    /**
     * ConfigurationError is thrown when invalid data is encountered while
     * parsing the configuration file.
     *
     * \ingroup grpconfig
     * \ingroup grpexceptions
     */
    class ConfigurationError :
        public Exception
    {
        public:
            /**
             * Constructor.
             *
             * \param line The configuration file's line that cause the exception.
             */
            ConfigurationError(const std::string & line);
    };

    /**
     * Configuration is used to obtain user configured settings.
     *
     * \ingroup grpconfig
     */
    class Configuration
    {
        private:
            struct Implementation;

            /// Our implementation.
            Implementation * _imp;

            /// \name Unwanted operations
            /// \{

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            Configuration(const Configuration & other);

            /// Unwanted assignment-operator: Do not implement. See EffCpp, Item 27.
            const Configuration & operator= (const Configuration & other);

            /// \}

            /**
             * Constructor.
             *
             * For internal use only.
             */
            Configuration();

            /// Read in the configuration file.
            void _read();

        public:
            /// Return the only instance of Configuration.
            static Configuration * instance();

            /**
             * \name Value Retrieval
             * \{
             * Return the value for a named configuration entry.
             *
             * \param name The name of the configuration entry that shall be
             * looked up.
             * \param default The default value that shall be returned if there
             * is no configuration entry that matchs the name.
             */

            int get_value(const std::string & name, int default_value);

            std::string get_value(const std::string & name, const std::string & default_value);

            /// \}

            /**
             * \name Value Setting
             * \{
             * Set the value for a named configuration entry.
             *
             * \param name The name of the configuration entry that shall be
             * modified.
             * \param value The value that shall be used.
             */

            void set_value(const std::string & name, int value);

            void set_value(const std::string & name, const std::string & value);

            /// \}

            /// Re-read the configuration file.
            void reread();
    };
}

#endif
