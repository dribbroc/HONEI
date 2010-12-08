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

#pragma once
#ifndef LIBUTIL_GUARD_STRING_TOKENIZER_HH
#define LIBUTIL_GUARD_STRING_TOKENIZER_HH 1

#include <honei/util/exception.hh>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

namespace honei
{
    std::vector<std::string> string_tokenizer(const std::string & str, const std::string & delimiters)
    {
        std::vector<std::string> tokens;
        std::string::size_type delimPos = 0, tokenPos = 0, pos = 0;

        if(str.length()<1)  return tokens;
        while(1)
        {
            delimPos = str.find_first_of(delimiters, pos);
            tokenPos = str.find_first_not_of(delimiters, pos);

            if(std::string::npos != delimPos)
            {
                if(std::string::npos != tokenPos)
                {
                    if(tokenPos<delimPos)
                    {
                        tokens.push_back(str.substr(pos,delimPos-pos));
                    }
                    else
                    {
                        //tokens.push_back("");
                    }
                }
                else
                {
                    //tokens.push_back("");
                }
                pos = delimPos+1;
            }
            else
            {
                if(std::string::npos != tokenPos)
                {
                    tokens.push_back(str.substr(pos));
                }
                else
                {
                    //tokens.push_back("");
                }
                break;
            }
        }

        return tokens;
    }
}
#endif
