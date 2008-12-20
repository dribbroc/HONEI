/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LBM_SCENARIO_CONTROLLER_BASE_HH
#define LBM_SCENARIO_CONTROLLER_BASE_HH

class ScenarioControllerBase
{
    public:
        virtual ~ScenarioControllerBase()
        {
        }

        static int get_precision(int scen_id);

        virtual void init(void) = 0;

        virtual void do_timestep(void) = 0;

        virtual void render(bool show_ground, bool use_quads, bool enable_alpha_blending, bool show_water, float alpha) = 0;
};
#endif
