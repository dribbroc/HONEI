/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based in parts upon 'exception.cc' from Paludis, which is:
 *     Copyright (c) 2005, 2006, 2007 Ciaran McCreesh
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

#include <honei/util/assertion.hh>
#include <honei/util/attributes.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/stringify.hh>

#include <iostream>

#include <cstdlib>
#include <cstring>
#include <cxxabi.h>
#include <list>
#include <string>

namespace
{
    std::string
    join(std::list<std::string>::const_iterator begin, std::list<std::string>::const_iterator end,
            const std::string & delimiter)
    {
        std::string result;

        if (begin != end)
            while (true)
            {
                result += *begin;
                if (++begin == end)
                    break;

                result += delimiter;
            }

        return result;
    }
}

namespace honei
{
    HONEI_THREAD_LOCAL std::list<std::string> * context_stack = 0;

    struct Exception::ContextData
    {
        std::list<std::string> local_context_stack;

        ContextData()
        {
            if (context_stack)
            {
                local_context_stack.assign(context_stack->begin(), context_stack->end());
            }
        }

        std::string backtrace(const std::string & delimiter) const
        {
            if (context_stack)
            {
                return join(local_context_stack.begin(), local_context_stack.end(), delimiter);
            }
            else return "";
        }
    };
}

using namespace honei;

Context::Context(const char * const file, const long line, const std::string & context)
{
    if (! context_stack)
    {
        context_stack = new std::list<std::string>;
    }

    context_stack->push_back(context + " (" + stringify(file) + ":" + stringify(line) +")");
}

Context::~Context()
{
    if (! context_stack)
        throw InternalError("no context!");

    context_stack->pop_back();

    if (context_stack->empty())
    {
        delete context_stack;
        context_stack = 0;
    }
}

std::string
Context::backtrace(const std::string & delimiter)
{
    if (! context_stack)
        return "";

    return join(context_stack->begin(), context_stack->end(), delimiter);
}


Exception::Exception(const std::string & message) throw () :
    _context_data(new ContextData),
    _message(message)
{
}

Exception::Exception(const Exception & other) :
    std::exception(other),
    _context_data(new ContextData(*other._context_data)),
    _message(other._message)
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
            _what_str += " (" + message() + ")";
            std::free(name);
        }
    }
    if (_what_str.empty())
        //_what_str = stringify(std::exception::what());
        _what_str += " (" + message() + ")";
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

PThreadError::PThreadError(const std::string & function, int errno) throw () :
    ExternalError("libpthread", function + " failed, " + stringify(strerror(errno)))
{
}

