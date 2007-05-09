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

#include "exception.hh"
#include "stringify.hh"

#include <cxxabi.h>

using namespace pg512;

Exception::Exception(const std::string & message) throw () :
    _message(message)
{
}

Exception::Exception(const Exception & other) :
    std::exception(other),
    _message(other._message)
{
}

Exception::~Exception() throw ()
{
}

const std::string &
Exception::message() const
{
    return _message;
}

const char *
Exception::what() const throw ()
{
#ifdef HAVE_CXA_DEMANGLE
    if (_what_str.empty())
    {
        int status(0);
        char * const name(abi::__cxa_demangle(
                    ("_Z" + stringify(std::exception::what())).c_str(), 0, 0, &status));
        if (0 == status)
        {
            _what_str = name;
            std::free(name);
        }
    }
#endif
    if (_what_str.empty())
        _what_str = stringify(std::exception::what());
    return _what_str.c_str();
}
