/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/libswe/relax_solver.hh>
#include <iostream>
#include <scenario_controller.hh>

int main(int argc, char ** argv)
{
#if defined (HONEI_SSE)
    ScenarioController<tags::CPU::SSE, float> controller(0);
#elif defined (HONEI_CELL)
    ScenarioController<tags::Cell, float> controller(0);
#else
    ScenarioController<tags::CPU, float> controller(0);
#endif

    controller.init();
    controller.do_timestep();
}
