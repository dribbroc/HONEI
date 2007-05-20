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

#include <libutil/log.hh>

#include <list>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

using namespace pg512;

class Log::LogOutput
{
    private:
        /// Our log level.
        LogLevel _level;

        /// Our output stream.
        std::ostream * _output;

    public:
        /**
         * Constructor.
         *
         * \param level Log level for this output.
         * \param name Filename for this output.
         **/
        LogOutput(LogLevel level, const std::string & name) :
            _level(level),
            _output(new std::fstream(name.c_str(), std::fstream::out))
        {
        }

        /**
         * Constructor.
         *
         * \param level Log level for this output.
         * \param name Stream for this output.
         **/
        LogOutput(LogLevel level, std::ostream * output) :
            _level(level),
            _output(output)
        {
        }

        /**
         * Out output operator.
         *
         * \param message Message to be logged.
         **/
        std::ostream & operator<< (const std::string & message)
        {
            return (*_output << message);
        }

        /// Returns our loglevel.
        const LogLevel level() const
        {
            return _level;
        }
};

Log::Log()
{
    _outputs.push_back(LogOutput(ll_minimal, &std::cerr));
}

Log *
Log::instance()
{
    static Log result;

    return &result;
}

void
Log::message(const LogLevel level, const std::string & message)
{
    for (std::list<LogOutput>::iterator o(_outputs.begin()), o_end(_outputs.end()) ;
            o != o_end ; ++o)
    {
        if (level <= o->level())
            (*o) << message << std::endl;
    }
}
