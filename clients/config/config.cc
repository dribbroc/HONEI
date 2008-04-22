/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of HONEI. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/configuration.hh>

#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace h = honei;

int main(int argc, char ** argv)
{
    std::cout << "# Configuration file '" << h::Configuration::instance()->filename() << "'" << std::endl;

    for (h::Configuration::ConstIterator i(h::Configuration::instance()->begin()),
            i_end(h::Configuration::instance()->end()) ; i != i_end ; ++i)
    {
        std::cout << std::setw(40) << i->first << std::setw(0) << " = "
            << i->second << std::endl;
    }

    return EXIT_SUCCESS;
}
