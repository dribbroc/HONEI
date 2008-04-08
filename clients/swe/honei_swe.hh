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

#ifndef SWE_HONEI_SWE_HH
#define SWE_HONEI_SWE_HH

#include <scenario_controller.hh>

ScenarioController<tags::CPU::SSE, float> * controller_f;
ScenarioController<tags::CPU::SSE, double> * controller_d;
void switch_scenario(int id);
void display();
void resize(int width, int height);
void keyboard(unsigned char key, int x, int y);
void keyboard_s(int key, int x, int y);
void mouse(int button, int state, int x, int y);
void menu_rendering(GLint index);
void menu_scenario(GLint index);
void menu_main(GLint index);
void ogl_init();


#endif
