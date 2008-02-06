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

#include <honei/libutil/configuration.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/mutex.hh>

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <sys/stat.h>

using namespace honei;

ConfigurationError::ConfigurationError(const std::string & line) :
    Exception("Invalid line in configuration file: '" + line + "'")
{
}

struct Configuration::Implementation
{
    typedef std::map<std::string, int> IntValueMap;

    typedef std::map<std::string, std::string> StringValueMap;

    /// Our mutex.
    Mutex * const mutex;

    /// Our map from names to integer values.
    IntValueMap int_value_map;

    /// Our map from names to string values.
    StringValueMap string_value_map;

    /// Our configuration file's name.
    std::string filename;

    /// Constructor.
    Implementation() :
        mutex(new Mutex)
    {
    }

    /// Destructor.
    ~Implementation()
    {
        delete mutex;
    }
};


Configuration::Configuration() :
    _imp(new Implementation)
{
    _read();
}

Configuration *
Configuration::instance()
{
    static Configuration result;

    return &result;
}

void
Configuration::_read()
{
    char * envvar(std::getenv("HONEI_CONFIG"));
    struct stat stat_info;

    if (! envvar)
    {
        if (0 != ::lstat("./honeirc", &stat_info))
        {
            envvar = std::getenv("HOME");

            if (envvar)
            {
                _imp->filename = std::string(envvar);
                _imp->filename += "/.honeirc";
            }
        }
        else
        {
            _imp->filename = "./honeirc";
        }
    }
    else
    {
        _imp->filename = envvar;
    }

    if (0 == ::lstat(_imp->filename.c_str(), &stat_info))
    {
        CONTEXT("When reading configuration file '" + _imp->filename + "'");

        std::fstream file(_imp->filename.c_str(), std::ios_base::in);
        std::string line;

        while (std::getline(file, line))
        {
            if (line.empty())
                continue;

            if ('#' == line.at(0))
                continue;

            std::string::size_type pos;
            if (std::string::npos == (pos = line.find_first_of("=")))
                throw ConfigurationError(line);

            std::string key(line.substr(0, pos));
            std::string::size_type last_non_whitespace(key.find_last_not_of(" \t"));
            if (key.size() - 1 != last_non_whitespace)
                key.erase(last_non_whitespace + 1);

            int int_value(0);
            std::string string_value("");
            if (pos + 1 < line.size())
            {
                string_value = line.substr(pos + 1);
                std::string::size_type first_non_whitespace(string_value.find_first_not_of(" \t"));
                string_value.erase(0, first_non_whitespace);

                std::stringstream value_stream(string_value);
                value_stream >> int_value;
            }

            _imp->int_value_map.insert(std::make_pair(key, int_value));
            _imp->string_value_map.insert(std::make_pair(key, string_value));
        }

        file.close();
    }
}

int
Configuration::get_value(const std::string & name, int default_value)
{
    Lock l(*_imp->mutex);
    Implementation::IntValueMap::const_iterator v(_imp->int_value_map.find(name));
    int result(default_value);

    if (v != _imp->int_value_map.end())
    {
        result = v->second;
    }

    return result;
}

std::string
Configuration::get_value(const std::string & name, const std::string & default_value)
{
    Lock l(*_imp->mutex);
    Implementation::StringValueMap::const_iterator v(_imp->string_value_map.find(name));
    std::string result(default_value);

    if (v != _imp->string_value_map.end())
    {
        result = v->second;
    }

    return result;
}

void
Configuration::set_value(const std::string & name, int value)
{
    Lock l(*_imp->mutex);
    Implementation::IntValueMap::iterator v(_imp->int_value_map.find(name));

    if (v != _imp->int_value_map.end())
    {
        v->second = value;
    }
    else
    {
        v = _imp->int_value_map.insert(Implementation::IntValueMap::value_type(name, value)).first;
    }
}

void
Configuration::set_value(const std::string & name, const std::string & value)
{
    Lock l(*_imp->mutex);
    Implementation::StringValueMap::iterator v(_imp->string_value_map.find(name));

    if (v != _imp->string_value_map.end())
    {
        v->second = value;
    }
    else
    {
        v = _imp->string_value_map.insert(Implementation::StringValueMap::value_type(name, value)).first;
    }
}

Configuration::ConstIterator
Configuration::begin() const
{
    Lock l(*_imp->mutex);
    const Implementation & imp(*_imp);

    return ConstIterator(imp.string_value_map.begin());
}

Configuration::ConstIterator
Configuration::end() const
{
    Lock l(*_imp->mutex);
    const Implementation & imp(*_imp);

    return ConstIterator(imp.string_value_map.end());
}

void
Configuration::reread()
{
    Lock l(*_imp->mutex);

    _imp->int_value_map.clear();
    _imp->string_value_map.clear();

    _read();
}

std::string
Configuration::filename() const
{
    Lock l(*_imp->mutex);

    return _imp->filename;
}
