/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#ifndef LBM_HONEI_POISSON_HH
#define LBM_HONEI_POISSON_HH

#include <scenario_controller_base.hh>
#include <honei/util/time_stamp.hh>

ScenarioControllerBase * controller_f;
ScenarioControllerBase * controller_d;

double rotation_x_increment;
double rotation_y_increment;
double rotation_z_increment;

double translation_x_increment;
double translation_y_increment;
double translation_z_increment;

bool calc;
bool filling;
bool use_quads;
bool show_ground;
bool show_water;
bool enable_shading;
bool enable_alpha_blending;
bool paused;
bool fullscreen;

float alpha;
int screen_width;
int screen_height;

double rotation_x;
double rotation_y;
double rotation_z;

double translation_x;
double translation_y;
double translation_z;

TimeStamp actual, last;

void switch_scenario(int id);
static void resize(int width, int height);
static void keyboard(unsigned char key, int x, int y);
static void keyboard_s(int key, int x, int y);
static void mouse(int button, int state, int x, int y);
static void menu_rendering(GLint index);
static void menu_scenario(GLint index);
static void menu_main(GLint index);
static void ogl_init();
static void display();
static void display_null(){}


#endif
