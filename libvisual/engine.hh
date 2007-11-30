/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LIBSWE_GUARD_ENGINE_HH
#define LIBSWE_GUARD_ENGINE_HH 1
/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <GL/glut.h>
#include <libla/dense_matrix.hh>
namespace honei
{
    namespace gl_globals {

        double rotation_x_increment = 0;
        double rotation_y_increment = 0;
        double rotation_z_increment = 0;

        double translation_x_increment = 0;
        double translation_y_increment = 0;
        double translation_z_increment = 0;

        bool filling = 0;

        int screen_width = 640;
        int screen_height = 480;

        double rotation_x = 0;
        double rotation_y = 0;
        double rotation_z = 0;

        double translation_x = 0;
        double translation_y = 0;
        double translation_z = 0;
    }

    class Engine
    {

        public:

            static void init(void)
            {
                glClearColor(0.0, 0.0, 0.2, 0.0);
                glShadeModel(GL_SMOOTH);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f,(GLfloat)gl_globals::screen_width/(GLfloat)gl_globals::screen_height,1.0f,1000.0f);
                glEnable(GL_DEPTH_TEST);
                glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
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
                    case 'r': case 'R':
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
                }
            }

            static void display(void)
            {
                int l_index;
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();

                glTranslatef(-2.0,-2.0,-20);
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


                //hardcoded test mesh: TODO: remove
                DenseMatrix<double> scalarfield(10 , 10 , double(1));
                glBegin(GL_QUADS);

                for(unsigned int i = 0; i < 9; ++i)
                {
                    for(unsigned int j = 0; j < 9; ++j)
                    {
                        glColor3f(0.0, 0.0, 1.0);
                        glVertex3d(i,j,scalarfield[i][j]);
                        glColor3f(0.0, 1.0, 1.0);
                        glVertex3d(i+1,j,scalarfield[i+1][j]);
                        glVertex3d(i+1,j+1,scalarfield[i+1][j+1]);
                        glVertex3d(i,j+1,scalarfield[i][j+1]);
                    }
                }
                glEnd();
                glFlush();
                glutSwapBuffers();

            }
    };

}
#endif
