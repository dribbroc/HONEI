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

#include <libutil/assertion.hh>
#include <libutil/exception.hh>
#include <libutil/stringify.hh>

#include <libebt/libebt.hh>

#include <cxxabi.h>

using namespace honei;

namespace
{
    struct ContextTag;
}

struct Context::ContextData
{
    libebt::BacktraceContext<ContextTag> context;

    ContextData(const char * const file, const long line, const std::string & s) :
        context(stringify(file) + ":" + stringify(line) + ": " + s)
    {
    }
};

Context::Context(const char * const file, const long line, const std::string & context) :
    _context_data(new ContextData(file, line, context))
{
}

Context::~Context()
{
    delete _context_data;
}

std::string
Context::backtrace(const std::string & delimiter)
{
    return libebt::BacktraceContext<ContextTag>::backtrace(delimiter);
}

struct Exception::ContextData :
    public libebt::Backtraceable<ContextTag>
{
};

Exception::Exception(const std::string & message) throw () :
    _message(message),
    _context_data(new ContextData)
{
}

Exception::Exception(const Exception & other) :
    std::exception(other),
    _message(other._message),
    _context_data(new ContextData(*other._context_data))
{
}

Exception::~Exception() throw ()
{
    delete _context_data;
}

const std::string &
Exception::message() const throw ()
{
    return _message;
}

std::string
Exception::backtrace(const std::string & delimiter) const
{
    return _context_data->backtrace(delimiter);
}

const char *
Exception::what() const throw ()
{
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
    if (_what_str.empty())
        _what_str = stringify(std::exception::what());
    return _what_str.c_str();
}

InternalError::InternalError(const std::string & message) throw () :
    Exception("Internal error: " + message)
{
}

ExternalError::ExternalError(const std::string & library, const std::string & message) throw () :
    Exception("External error in '" + library + "': " + message)
{
}
