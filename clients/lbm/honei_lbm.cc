/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <GL/glut.h>
#include <iostream>
#include <scenario_controller_base.hh>
#include <scenario_controller.hh>
#include <scenario_controller_grid.hh>
#include <scenario_controller_dat.hh>
#include <honei_lbm.hh>

int main(int argc, char ** argv)
{
    //OGL
    rotation_x_increment = 0;
    rotation_y_increment = 0;
    rotation_z_increment = 0;

    translation_x_increment = 0;
    translation_y_increment = 0;
    translation_z_increment = 0;

    calc = true;
    filling = 1;
    use_quads = true;
    show_ground = true;
    show_water = true;
    enable_shading = true;
    enable_alpha_blending = true;
    paused = false;
    fullscreen = false;

    alpha = 0.8;
    screen_width = 640;
    screen_height = 480;

    rotation_x = -25;
    rotation_y = 0;
    rotation_z = 0;

    translation_x = 0;
    translation_y = 0;
    translation_z = 0;

    int i =1;
    int * pi = &i;

    char * c = "Visual LBM";
    char ** cp = &c;
    glutInit(pi,cp);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
    glutInitWindowSize(screen_width, screen_height);
    glutInitWindowPosition(0,0);
    glutCreateWindow("HONEI SWE");
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutReshapeFunc(resize);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(keyboard_s);
    glutMouseFunc(mouse);
    GLint menu_id_rendering = glutCreateMenu(menu_rendering);
    glutAddMenuEntry("Toggle fill mode", 2);
    glutAddMenuEntry("Toggle ground", 3);
    glutAddMenuEntry("Toggle water", 4);
    glutAddMenuEntry("Toggle shading", 5);
    glutAddMenuEntry("Toggle primitive type", 6);
    glutAddMenuEntry("Toggle alpha blending", 7);
    GLint menu_id_scenario_mono = glutCreateMenu(menu_scenario);
    glutAddMenuEntry("Laminar flow: Circular dam break 50x50", 0);
    glutAddMenuEntry("Laminar flow: Circular dam break above uneven bed 50x50", 1);
    GLint menu_id_scenario_grid = glutCreateMenu(menu_scenario);
    glutAddMenuEntry("Laminar flow: Circular dam break 50x50 ", 100);
    glutAddMenuEntry("Laminar flow: Circular dam break 50x50 with cuboidal obstacles", 101);
    glutAddMenuEntry("Laminar flow: Partial dam break 50x50", 102);
    glutAddMenuEntry("Laminar flow: Circular dam break above uneven bed 50x50", 103);
    glutAddMenuEntry("Laminar flow: Circular dam break above uneven bed (b) 100x200", 104);
    GLint menu_id_scenario_multi = glutCreateMenu(menu_scenario);
    glutAddMenuEntry("Laminar flow: Circular dam break 50x50 ", 1100);
    glutAddMenuEntry("Laminar flow: Circular dam break 50x50 with cuboidal obstacles", 1101);
    glutAddMenuEntry("Laminar flow: Partial dam break 50x50", 1102);
    glutAddMenuEntry("Laminar flow: Circular dam break above uneven bed 50x50", 1103);
    glutAddMenuEntry("Laminar flow: Circular dam break above uneven bed (b) 100x200", 1104);
    GLint menu_id_scenario = glutCreateMenu(menu_scenario);
    glutAddSubMenu("Mono", menu_id_scenario_mono);
    glutAddSubMenu("Grid", menu_id_scenario_grid);
    glutAddSubMenu("Multi", menu_id_scenario_multi);
    glutAddMenuEntry("Animate GnuPlot dat files", -1);
    GLint menu_id_main = glutCreateMenu(menu_main);
    glutAddMenuEntry("Restart scenario", 0);
    glutAddSubMenu("Rendering", menu_id_rendering);
    glutAddSubMenu("Scenarios", menu_id_scenario);
    glutAddMenuEntry("Exit programm", 10);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    glutMotionFunc(NULL);
    glutPassiveMotionFunc(NULL);
    ogl_init();


    controller_f = 0;
    controller_d = 0;

    switch_scenario(100);

    glutMainLoop();


    delete controller_f;
    delete controller_d;
}

void switch_scenario(int id)
{
    std::cout<<"switching to scenario: "<<id<<std::endl;
    glutIdleFunc(NULL);
    glutDisplayFunc(display_null);
    calc = false;

    if (id == -1)
    {
            delete controller_f;
            delete controller_d;
            controller_f = 0;
            controller_d = new ScenarioControllerDat<tags::CPU, double> (id);
            controller_d->init();
    }
    else if (id < 100)
    {
        if (ScenarioController<tags::CPU, float>::get_precision(id) == 0)
        {
            delete controller_f;
            delete controller_d;
            controller_d = 0;
            controller_f = new ScenarioController<tags::CPU, float> (id);
            controller_f->init();
        }
        else if (ScenarioController<tags::CPU, float>::get_precision(id) == 1)
        {
            delete controller_f;
            delete controller_d;
            controller_f = 0;
            controller_d = new ScenarioController<tags::CPU, double> (id);
            controller_d->init();
        }
    }
    else if (id < 200)
    {
        if (ScenarioControllerGrid<tags::CPU, float>::get_precision(id) == 0)
        {
            delete controller_f;
            delete controller_d;
            controller_d = 0;
            controller_f = new ScenarioControllerGrid<tags::CPU, float> (id);
            controller_f->init();
        }
        else if (ScenarioControllerGrid<tags::CPU, float>::get_precision(id) == 1)
        {
            delete controller_f;
            delete controller_d;
            controller_f = 0;
            controller_d = new ScenarioControllerGrid<tags::CPU, double> (id);
            controller_d->init();
        }
    }
    else if (id < 1200)
    {
        if (ScenarioControllerGrid<tags::CPU, float>::get_precision(id) == 0)
        {
            delete controller_f;
            delete controller_d;
            controller_d = 0;
            controller_f = new ScenarioControllerGrid<tags::CPU::MultiCore, float> (id - 1000);
            controller_f->init();
        }
        else if (ScenarioControllerGrid<tags::CPU, float>::get_precision(id) == 1)
        {
            delete controller_f;
            delete controller_d;
            controller_f = 0;
            controller_d = new ScenarioControllerGrid<tags::CPU::MultiCore, double> (id - 1000);
            controller_d->init();
        }
    }
    calc = true;
    glutDisplayFunc(display);
    glutIdleFunc(display);
}

static void resize (int width, int height)
{
    screen_width=width;
    screen_height=height;
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0,0,screen_width, screen_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,1000.0f);
    glutPostRedisplay ();
}

static void keyboard (unsigned char key, int x, int y)
{
    switch (key)
    {
        case ' ':
            rotation_x_increment=0;
            rotation_y_increment=0;
            rotation_z_increment=0;
            translation_x_increment=0;
            translation_y_increment=0;
            translation_z_increment=0;

            break;
        case 'r':
            if (filling==0)
            {
                glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
                filling=1;
            }
            else
            {
                glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
                filling=0;
            }
            break;
        case 'g':
            show_ground = !show_ground;
            break;
        case 'w':
            show_water = !show_water;
            break;
        case 's':
            enable_shading = !enable_shading;
            break;
        case 'p':
            if (paused)
            {
                paused = false;
                glutIdleFunc(display);
            }
            else
            {
                paused = true;
                glutIdleFunc(NULL);
            }
            break;
        case 'n':
            translation_x_increment = translation_x_increment -0.1;
            break;

        case 'm':
            translation_x_increment = translation_x_increment +0.1;
            break;

        case ',':
            translation_y_increment = translation_y_increment -0.1;
            break;

        case '.':
            translation_y_increment = translation_y_increment +0.1;
            break;

        case 'v':
            translation_z_increment = translation_z_increment -0.1;
            break;

        case 'b':
            translation_z_increment = translation_z_increment +0.1;
            break;

        case 'q':
            exit (0);
            break;

        case char(27):
            exit(0);
            break;
    }
}

static void keyboard_s (int key, int x, int y)
{
    switch (key)
    {
        case GLUT_KEY_UP:
            rotation_x_increment = rotation_x_increment +0.5;
            break;
        case GLUT_KEY_DOWN:
            rotation_x_increment = rotation_x_increment -0.5;
            break;
        case GLUT_KEY_LEFT:
            rotation_y_increment = rotation_y_increment +0.5;
            break;
        case GLUT_KEY_RIGHT:
            rotation_y_increment = rotation_y_increment -0.5;
            break;
        case GLUT_KEY_F5:
            ogl_init();
            if (controller_f)
            {
                controller_f->init();
            }
            else if (controller_d)
            {
                controller_d->init();
            }
            break;
        case GLUT_KEY_F6:
            rotation_x = -25;
            rotation_y = 0;
            rotation_z = 0;

            translation_x = 0;
            translation_y = 0;
            translation_z = 0;
            break;
        case GLUT_KEY_F8:
            if (fullscreen)
            {
                glutReshapeWindow(640,480);
                glutPositionWindow(0,0);
                glutPostRedisplay();
                fullscreen = false;
            }
            else
            {
                glutFullScreen();
                fullscreen = true;
            }
            break;
    }
}

static void mouse (int button, int state, int x, int y)
{
    /*if (button == GLUT_LEFT_BUTTON && state ==GLUT_DOWN)
      glutSetCursor(GLUT_CURSOR_INFO);
      if (button == GLUT_LEFT_BUTTON && state ==GLUT_UP)
      glutSetCursor(GLUT_CURSOR_INHERIT);*/
}

static void menu_rendering(GLint index)
{
    switch (index)
    {
        case 2:
            if (filling==0)
            {
                glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
                filling=1;
            }
            else
            {
                glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
                filling=0;
            }
            break;
        case 3:
            show_ground = !show_ground;
            break;
        case 4:
            show_water = !show_water;
            break;
        case 5:
            enable_shading = !enable_shading;
            break;
        case 6:
            use_quads = !use_quads;
            break;
        case 7:
            enable_alpha_blending = !enable_alpha_blending;
            break;
    }
}

static void menu_scenario(GLint index)
{
    switch_scenario(index);
}

static void menu_main(GLint index)
{
    switch(index)
    {
        case 0:
            ogl_init();
            if (controller_f)
            {
                controller_f->init();
            }
            else if (controller_d)
            {
                controller_d->init();
            }
            break;
        case 10:
            exit(0);
            break;
    }
}

static void ogl_init(void)
{
    last.take();
    glClearColor(0.0, 0.0, 0.2, 0.0);
    if (enable_shading) glShadeModel(GL_SMOOTH);
    else glShadeModel(GL_FLAT);
    glViewport(0,0,screen_width, screen_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,1000.0f);
    glEnable(GL_DEPTH_TEST);
    if (filling) glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    else glPolygonMode (GL_FRONT_AND_BACK, GL_LINES);
    //eye candy
    //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    //glEnable(GL_POLYGON_SMOOTH);
}

static void display(void)
{
    if (enable_shading) glShadeModel(GL_SMOOTH);
    else glShadeModel(GL_FLAT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef(-2.0,-10.0,-45);
    rotation_x = rotation_x + rotation_x_increment;
    rotation_y = rotation_y + rotation_y_increment;
    rotation_z = rotation_z + rotation_z_increment;
    if (rotation_x > 359) rotation_x = 0;
    if (rotation_y > 359) rotation_y = 0;
    if (rotation_z > 359) rotation_z = 0;
    glRotatef(rotation_x,1.0,0.0,0.0);
    glRotatef(rotation_y,0.0,1.0,0.0);
    glRotatef(rotation_z,0.0,0.0,1.0);
    //new:
    glRotatef(-45.0f,1.0, 0.0, 0.0);
    glRotatef(45.0f,0.0, 0.0, 1.0);
    translation_x = translation_x + translation_x_increment;
    translation_y = translation_y + translation_y_increment;
    translation_z = translation_z + translation_z_increment;

    glTranslatef(translation_x, 0.0, 0.0);
    glTranslatef(0.0, translation_y, 0.0);
    glTranslatef(0.0, 0.0 , translation_z);

    do
    {
        actual.take();
    }
    while(actual.usec() - last.usec() < 60000ul); // 1/25 = 40000
    last.take();

    if (calc)
    {
        if (controller_f)
        {
            controller_f->do_timestep();
            controller_f->render(show_ground, use_quads, enable_alpha_blending, show_water, alpha);
        }
        else if (controller_d)
        {
            controller_d->do_timestep();
            controller_d->render(show_ground, use_quads, enable_alpha_blending, show_water, alpha);
        }
    }

    glutSwapBuffers();

}
