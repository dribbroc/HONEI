/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_EXCEPTION_HH
#define LIBUTIL_GUARD_EXCEPTION_HH 1

#include <string>

namespace pg512 ///< \todo Namespace name?
{
    class Exception :
        public std::exception
    {
        private:
            const std::string _message;

        protected:
            Exception(const std::string & message) throw ();

            Exception(const Exception & e);

        public:
            virtual ~Exception() throw ();

            const std::string & message() const;
    };
}

#endif
