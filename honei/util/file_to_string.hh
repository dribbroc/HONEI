
/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_FILE_TO_STRING_HH
#define LIBUTIL_GUARD_FILE_TO_STRING_HH 1

#include <honei/util/exception.hh>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
namespace honei
{
    std::string file_to_string(std::string file)
    {
        size_t size;
        char*  str;
        std::string s;
        const char * filename = file.c_str();

        std::fstream f(filename, (std::fstream::in | std::fstream::binary));

        if(f.is_open())
        {
            size_t fileSize;
            f.seekg(0, std::fstream::end);
            size = fileSize = f.tellg();
            f.seekg(0, std::fstream::beg);

            str = new char[size+1];
            if(!str)
            {
                f.close();
                return NULL;
            }

            f.read(str, fileSize);
            f.close();
            str[size] = '\0';

            s = str;

            return s;
        }
        else throw InternalError("file_to_string: Could not open " + file);
    }
}
#endif
