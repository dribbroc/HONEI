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
#include <honei_swe.hh>
#include <scenario_controller.hh>

int main(int argc, char ** argv)
{
    controller_f = 0;
    controller_d = 0;

    switch_scenario(0);

    if (controller_f)
    {
        controller_f->init();
        controller_f->do_timestep();
    }
    else if (controller_d)
    {
        controller_d->init();
        controller_d->do_timestep();
    }



    delete controller_f;
    delete controller_d;
}

void switch_scenario(int id)
{
#if defined (HONEI_SSE)
    if (ScenarioController<tags::CPU, float>::get_precision(id) == 0)
    {
        delete controller_f;
        delete controller_d;
        controller_d = 0;
        controller_f = new ScenarioController<tags::CPU::SSE, float> (id);
    }
    else if (ScenarioController<tags::CPU, float>::get_precision(id) == 1)
    {
        delete controller_f;
        delete controller_d;
        controller_f = 0;
        controller_d = new ScenarioController<tags::CPU::SSE, double> (id);
    }
#endif
}
