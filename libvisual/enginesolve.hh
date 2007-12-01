/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LibVisual C++ library. LibVisual is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibVisual is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBSWE_GUARD_ENGINESOLVE_HH
#define LIBSWE_GUARD_ENGINESOLVE_HH 1

#include <GL/glut.h>
#include <libla/dense_matrix.hh>
#include <libswe/solver.hh>

namespace honei
{
    namespace gl_globals {

        double rotation_x_increment = 0;
        double rotation_y_increment = 0;
        double rotation_z_increment = 0;

        double translation_x_increment = 0;
        double translation_y_increment = 0;
        double translation_z_increment = 0;

        bool filling = 1;
        bool use_quads = true;
        bool show_ground = true;
        bool show_water = true;
        bool enable_shading = true;
        bool enable_alpha_blending = true;
        bool pause = false;
        bool fullscreen = false;

        float alpha = 0.8;
        int screen_width = 800;
        int screen_height = 600;

        double rotation_x = 0;
        double rotation_y = 0;
        double rotation_z = 0;

        double translation_x = 0;
        double translation_y = 0;
        double translation_z = 0;

        GLint menu_id_main;
        GLint menu_id_scenario;
        GLint menu_id_rendering;

        //globally defined solver:
        ulint dwidth =40;
        ulint dheight =40;
        //DenseMatrix<float> height(dheight, dwidth, float(5));
        DenseMatrix<float> height(dheight, dwidth, float(5));

        DenseMatrix<float> bottom(dheight, dwidth, float(0));
        DenseMatrix<float> u1(dheight, dwidth, float(4));
        DenseMatrix<float> u2(dheight, dwidth, float(4));
        unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
        DenseVector<float> u(entries, float(1));
        DenseVector<float> v(entries, float(1));
        DenseVector<float> w(entries, float(1));
        DenseVector<float> bx (entries/3, float(0));
        DenseVector<float> by (entries/3, float(0));
        DenseVector<float> c (3,float(5));
        DenseVector<float> d (3,float(5));
        float deltax = 5;
        float deltay = 5;
        float deltat = 5./24.;
        double eps = 10e-6;
        float manning = float(0);
        RelaxSolver<tags::CPU::SSE, float, float, float, float, float> solver( &height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);

    }

    class Engine
    {

        public:


            static void init(void)
            {
                glClearColor(0.0, 0.0, 0.2, 0.0);
                if (gl_globals::enable_shading) glShadeModel(GL_SMOOTH);
                else glShadeModel(GL_FLAT);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f,(GLfloat)gl_globals::screen_width/(GLfloat)gl_globals::screen_height,1.0f,1000.0f);
                glEnable(GL_DEPTH_TEST);
                if (gl_globals::filling) glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
                else glPolygonMode (GL_FRONT_AND_BACK, GL_LINES);
                //eye candy
                //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
                //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
                //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
                //glEnable(GL_POLYGON_SMOOTH);



                //Solver setup (Scenario):
                for(unsigned int i = 0; i < gl_globals::dheight; i++)
                {
                    if(i == gl_globals::dheight/2 -4 || i == gl_globals::dheight/2 +5)
                    {
                        gl_globals::height[i][gl_globals::dheight/2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +1] = float(15);
                    }
                    else if(i == gl_globals::dheight/2 -3 || i == gl_globals::dheight/2 +4)
                    {
                        gl_globals::height[i][gl_globals::dheight/2 -2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 ] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +3] = float(15);
                    }
                    else if(i == gl_globals::dheight/2 -2 || i == gl_globals::dheight/2 -1 || i == gl_globals::dheight/2 +2 || i == gl_globals::dheight/2 + 3)
                    {
                        gl_globals::height[i][gl_globals::dheight/2 -3] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +3] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +4] = float(15);

                    }
                    else if (i == gl_globals::dheight/2 || i == gl_globals::dheight/2 +1)
                    {
                        gl_globals::height[i][gl_globals::dheight/2 -4] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -3] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 -1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +1] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +2] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +3] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +4] = float(15);
                        gl_globals::height[i][gl_globals::dheight/2 +5] = float(15);

                    }
                }
                /*
                for(ulint i = 0; i< gl_globals::height.rows(); ++i)
                  {
                  for(ulint j=gl_globals::height.columns()-10; j<gl_globals::height.columns(); ++j)
                  {
                  gl_globals::height[i][j] = float(10);
                  }
                  }
                  for(ulint i = 0; i< gl_globals::bottom.rows(); ++i)
                  {
                  for(ulint j=0; j<gl_globals::bottom.columns()-10; ++j)
                  {
                  gl_globals::bottom[i][j] = float(1);
                  if(j>4 && j< gl_globals::bottom.columns()-9)
                  {
                  if(i <8 || i > 13)
                  gl_globals::bottom[i][j] = float(3);
                  else
                  gl_globals::bottom[i][j] = float(1);
                  }
                  }
                  }*//*
                  gl_globals::c[0] = 10;
                  gl_globals::c[1] = 6;
                  gl_globals::c[2] = 11;
                  gl_globals::d[0] = 10;
                  gl_globals::d[1] = 5;
                  gl_globals::d[2] = 11;*/

                  gl_globals::c[0] = 12;
                  gl_globals::c[1] = 7;
                  gl_globals::c[2] = 12;
                  gl_globals::d[0] = 12;
                  gl_globals::d[1] = 7;
                  gl_globals::d[2] = 12;

                  gl_globals::solver.do_preprocessing();
            }

            static void resize (int width, int height)
            {
                gl_globals::screen_width=width;
                gl_globals::screen_height=height;
                glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f,(GLfloat)gl_globals::screen_width/(GLfloat)gl_globals::screen_height,1.0f,1000.0f);
                glutPostRedisplay ();
            }

            static void mouse (int button, int state, int x, int y)
            {
                if (button == GLUT_LEFT_BUTTON && state ==GLUT_DOWN)
                    glutSetCursor(GLUT_CURSOR_INFO);
                if (button == GLUT_LEFT_BUTTON && state ==GLUT_UP)
                    glutSetCursor(GLUT_CURSOR_INHERIT);
            }

            static void keyboard (unsigned char key, int x, int y)
            {
                switch (key)
                {
                    case ' ':
                        gl_globals::rotation_x_increment=0;
                        gl_globals::rotation_y_increment=0;
                        gl_globals::rotation_z_increment=0;
                        gl_globals::translation_x_increment=0;
                        gl_globals::translation_y_increment=0;
                        gl_globals::translation_z_increment=0;

                        break;
                    case 'r':
                        if (gl_globals::filling==0)
                        {
                            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
                            gl_globals::filling=1;
                        }
                        else
                        {
                            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
                            gl_globals::filling=0;
                        }
                        break;
                    case 'g':
                        gl_globals::show_ground = !gl_globals::show_ground;
                        break;
                    case 'w':
                        gl_globals::show_water = !gl_globals::show_water;
                        break;
                    case 's':
                        gl_globals::enable_shading = !gl_globals::enable_shading;
                        break;
                    case 'p':
                        if (gl_globals::pause)
                        {
                            gl_globals::pause = false;
                            glutIdleFunc(Engine::display);
                        }
                        else
                        {
                            gl_globals::pause = true;
                            glutIdleFunc(NULL);
                        }
                        break;
                    case 'n':
                        gl_globals::translation_x_increment = gl_globals::translation_x_increment -0.1;
                        break;

                    case 'm':
                        gl_globals::translation_x_increment = gl_globals::translation_x_increment +0.1;
                        break;

                    case ',':
                        gl_globals::translation_y_increment = gl_globals::translation_x_increment -0.1;
                        break;

                    case '.':
                        gl_globals::translation_y_increment = gl_globals::translation_x_increment +0.1;
                        break;

                    case 'v':
                        gl_globals::translation_y_increment = gl_globals::translation_x_increment -0.1;
                        break;

                    case 'b':
                        gl_globals::translation_y_increment = gl_globals::translation_x_increment +0.1;
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
                        gl_globals::rotation_x_increment = gl_globals::rotation_x_increment +0.5;
                        break;
                    case GLUT_KEY_DOWN:
                        gl_globals::rotation_x_increment = gl_globals::rotation_x_increment -0.5;
                        break;
                    case GLUT_KEY_LEFT:
                        gl_globals::rotation_y_increment = gl_globals::rotation_y_increment +0.5;
                        break;
                    case GLUT_KEY_RIGHT:
                        gl_globals::rotation_y_increment = gl_globals::rotation_y_increment -0.5;
                        break;
                    case GLUT_KEY_F5:
                        init();
                        break;
                    case GLUT_KEY_F8:
                        if (gl_globals::fullscreen)
                        {
                            glutReshapeWindow(640,480);
                            glutPositionWindow(0,0);
                            glutPostRedisplay();
                            gl_globals::fullscreen = false;
                        }
                        else
                        {
                            glutFullScreen();
                            gl_globals::fullscreen = true;
                        }
                        break;
                }
            }

            static void menu_rendering(GLint index)
            {
                switch (index)
                {
                    case 2:
                        if (gl_globals::filling==0)
                        {
                            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
                            gl_globals::filling=1;
                        }
                        else
                        {
                            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
                            gl_globals::filling=0;
                        }
                        break;
                    case 3:
                        gl_globals::show_ground = !gl_globals::show_ground;
                        break;
                    case 4:
                        gl_globals::show_water = !gl_globals::show_water;
                        break;
                    case 5:
                        gl_globals::enable_shading = !gl_globals::enable_shading;
                        break;
                    case 6:
                        gl_globals::use_quads = !gl_globals::use_quads;
                        break;
                }
            }

            static void menu_scenario(GLint index)
            {
            }

            static void menu_main(GLint index)
            {
                switch(index)
                {
                    case 0:
                        init();
                        break;
                    case 10:
                        exit(0);
                        break;
                }
            }

            static void display(void)
            {
                //int l_index;
                if (gl_globals::enable_shading) glShadeModel(GL_SMOOTH);
                else glShadeModel(GL_FLAT);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();

                glTranslatef(-2.0,-10.0,-45);
                gl_globals::rotation_x = gl_globals::rotation_x + gl_globals::rotation_x_increment;
                gl_globals::rotation_y = gl_globals::rotation_y + gl_globals::rotation_y_increment;
                gl_globals::rotation_z = gl_globals::rotation_z + gl_globals::rotation_z_increment;
                if (gl_globals::rotation_x > 359) gl_globals::rotation_x = 0;
                if (gl_globals::rotation_y > 359) gl_globals::rotation_y = 0;
                if (gl_globals::rotation_z > 359) gl_globals::rotation_z = 0;
                glRotatef(gl_globals::rotation_x,1.0,0.0,0.0);
                glRotatef(gl_globals::rotation_y,0.0,1.0,0.0);
                glRotatef(gl_globals::rotation_z,0.0,0.0,1.0);
                //new:
                glRotatef(-45.0f,1.0, 0.0, 0.0);
                glRotatef(45.0f,0.0, 0.0, 1.0);
                gl_globals::translation_x = gl_globals::translation_x + gl_globals::translation_x_increment;
                gl_globals::translation_y = gl_globals::translation_y + gl_globals::translation_y_increment;
                gl_globals::translation_z = gl_globals::translation_z + gl_globals::translation_z_increment;

                glTranslatef(gl_globals::translation_x, 0.0, 0.0);
                glTranslatef(0.0, gl_globals::translation_y, 0.0);
                glTranslatef(0.0, 0.0 , gl_globals::translation_z);


                gl_globals::solver.solve();


                if (gl_globals::show_ground)
                {
                    if(gl_globals::use_quads)
                    {
                        glBegin(GL_QUADS);
                        for(unsigned int i = 0; i < gl_globals::dwidth-1; ++i)
                        {
                            for(unsigned int j = 0; j < gl_globals::dheight-1; ++j)
                            {
                                glColor3f(1.0, 0.0, 0.0);
                                glVertex3d(i,j,gl_globals::bottom[i][j]);
                                glColor3f(1.0, 0.8, 0.0);
                                glVertex3d(i+1,j,gl_globals::bottom[i+1][j]);
                                glVertex3d(i+1,j+1,gl_globals::bottom[i+1][j+1]);
                                glVertex3d(i,j+1,gl_globals::bottom[i][j+1]);
                            }
                        }
                        glEnd();
                    }
                    else
                    {
                        glBegin(GL_TRIANGLE_STRIP);
                        for(unsigned int i = 0; i < gl_globals::dwidth-1; ++i)
                        {
                            for(unsigned int j = 0; j < gl_globals::dheight; j++)
                            {
                                glColor3f(1.0, 0.0, 0.0);
                                glVertex3d(i,j,gl_globals::bottom[i][j]);
                                glColor3f(1.0, 0.8, 0.0);
                                glVertex3d(i+1,j,gl_globals::bottom[i+1][j]);
                            }
                            ++i;
                            if (i >= gl_globals::dwidth-1)
                                break;
                            for(int j2 = gl_globals::dheight-2; j2 >= 0; --j2)
                            {
                                glVertex3d(i,j2,gl_globals::bottom[i][j2]);
                                glColor3f(1.0, 0.0, 0.0);
                                glVertex3d(i+1,j2,gl_globals::bottom[i+1][j2]);
                            }
                        }
                        glEnd();
                    }
                }
                if(gl_globals::enable_alpha_blending)
                {
                    glEnable (GL_BLEND);
                    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                }

                if (gl_globals::show_water)
                {
                    if(gl_globals::use_quads)
                    {
                        glBegin(GL_QUADS);
                        for(unsigned int i = 0; i < gl_globals::dwidth-1; ++i)
                        {
                            for(unsigned int j = 0; j < gl_globals::dheight-1; ++j)
                            {
                                glColor4f(0.0, 0.0, 1.0, gl_globals::alpha);
                                glVertex3d(i,j,gl_globals::height[i][j] + gl_globals::bottom[i][j]);
                                glColor4f(0.0, 1.0, 1.0, gl_globals::alpha);
                                glVertex3d(i+1,j,gl_globals::height[i+1][j] + gl_globals::bottom[i+1][j]);
                                glVertex3d(i+1,j+1,gl_globals::height[i+1][j+1] + gl_globals::bottom[i+1][j+1]);
                                glVertex3d(i,j+1,gl_globals::height[i][j+1] + gl_globals::bottom[i][j+1]);
                            }
                        }
                        glEnd();
                    }
                    else
                    {
                        glBegin(GL_TRIANGLE_STRIP);
                        for(unsigned int i = 0; i < gl_globals::dwidth-1; ++i)
                        {
                            for(unsigned int j = 0; j < gl_globals::dheight; j++)
                            {
                                glColor3f(0.0, 1.0, 1.0);
                                glVertex3d(i,j,gl_globals::height[i][j] + gl_globals::bottom[i][j]);
                                glColor3f(0.0, 0.0, 1.0);
                                glVertex3d(i+1,j,gl_globals::height[i+1][j] + gl_globals::bottom[i+1][j]);
                            }
                            ++i;
                            if (i >= gl_globals::dwidth-1)
                                break;
                            for(int j2 = gl_globals::dheight-2; j2 >= 0; --j2)
                            {
                                glVertex3d(i,j2,gl_globals::height[i][j2] + gl_globals::bottom[i][j2]);
                                glColor3f(0.0, 1.0, 1.0);
                                glVertex3d(i+1,j2,gl_globals::height[i+1][j2] + gl_globals::bottom[i+1][j2]);
                            }
                        }
                        glEnd();
                    }
                }

                glutSwapBuffers();

            }
    };

}
#endif
