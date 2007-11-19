/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2005, 2006, 2007 Ciaran McCreesh <ciaranm@ciaranm.org>
 *
 * This file is part of the Paludis package manager. Paludis is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * Paludis is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBUTIL_GUARD_STRINGIFY_HH
#define LIBUTIL_GUARD_STRINGIFY_HH 1

#ifdef HONEI_CELL
# include <cell/cell.hh>
#endif

#include <sstream>
#include <string>
#include <tr1/memory>

namespace honei
{
    /**
     * For use by stringify.
     *
     * \ingroup grpstringify
     */
    namespace stringify_internals
    {
        /**
         * Check that T_ is a sane type to be stringified.
         *
         * \ingroup grpstringify
         */
        template <typename T_>
        struct CheckType
        {
            /// Yes, we are a sane type.
            enum { value = 0 } Value;
        };

        /**
         * Check that T_ is a sane type to be stringified.
         *
         * \ingroup grpstringify
         */
        template <typename T_>
        struct CheckType<T_ *>
        {
            /// Yes, we are a sane type.
            enum { value = 0 } Value;
        };

        /**
         * Check that T_ is a sane type to be stringified, which it isn't
         * if it's a shared_ptr.
         *
         * \ingroup grpstringify
         */
        template <typename T_>
        struct CheckType<std::tr1::shared_ptr<T_> >
        {
        };
    }

    /**
     * Convert item to a string.
     *
     * \ingroup grpstringify
     */
    template <typename T_>
    std::string
    stringify(const T_ & item)
    {
        /* check that we're not trying to stringify a pointer or somesuch */
        int check_for_stringifying_silly_things = stringify_internals::CheckType<T_>::value;

        std::ostringstream s;
        s << item;
        return s.str();
    }

    /**
     * Convert item to a string (overload for std::string).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const std::string & item)
    {
        return item;
    }

    /**
     * Convert item to a string (overload for char).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for unsigned char).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const unsigned char & item)
    {
        return std::string(1, item);
    }

    /**
     * Convert item to a string (overload for bool).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const bool & item)
    {
        return item ? "true" : "false";
    }

    /**
     * Convert item to a string (overload for char *, which isn't a
     * screwup like other pointers).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const char * const item)
    {
        return std::string(item);
    }

#ifdef HONEI_CELL
    /**
     * Convert item to a string (overload for EffectiveAddress).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const cell::EffectiveAddress & item)
    {
        std::ostringstream s;
        s << std::hex << item;
        return s.str();
    }

    /**
     * Convert item to a string (overload for LocalStoreAddress).
     *
     * \ingroup grpstringify
     */
    inline std::string
    stringify(const cell::LocalStoreAddress & item)
    {
        std::ostringstream s;
        s << std::hex << item;
        return s.str();
    }
#endif
}

#endif
