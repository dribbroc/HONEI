/* vim: set sw=4 sts=4 et nofoldenable nu : */

#ifndef LIBSWE_GUARD_ENGINE_HH
#define LIBSWE_GUARD_ENGINE_HH 1
 /*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>

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
#include <honei/libla/dense_matrix.hh>
#include <honei/graph/evolving_graph.hh>
#include <honei/graph/evolving_animator.hh>
#include <iostream>

namespace honei
{
    class Color
    {
        public:
            float r, g, b;

        Color(float r, float g, float b)
        {
            this->r = r;
            this->g = g;
            this->b = b;
        }
    };

    namespace gl_globals {
        typedef double DataType_;
        typedef tags::CPU Tag_;

        double rotation_x_increment = 0;
        double rotation_y_increment = 0;
        double rotation_z_increment = 0;

        double translation_x_increment = 0;
        double translation_y_increment = 0;
        double translation_z_increment = 0;

        bool filling = 0;

        int screen_width = 800;
        int screen_height = 600;

        double rotation_x = 0;
        double rotation_y = 0;
        double rotation_z = 0;

        double translation_x = 0;
        double translation_y = 0;
        double translation_z = 0;

        Color * colors[10];
        float ebene_z = 0.0f;
        bool calculate = false;
        void * animator;
        
        
        float edgeMaterial[] = {0.2, 0.2, 0.2, 0.8};
        float edgeTransparent[] = {0.4, 0.4, 0.6, 0.4};
    }
    
    template <typename Tag_, typename DataType_>
    class EngineEvolving
    {
        private:
        
        public:
            static void setTestCase(EvolvingAnimator<Tag_, DataType_> & animator)
            {               
                gl_globals::animator = &animator;
            }            
            
            static void init(void)
            {
                    gl_globals::colors[0] = new Color(1.0f, 0.0f, 0.0f);
                    gl_globals::colors[1] = new Color(0.0f, 1.0f, 0.0f);
                    gl_globals::colors[2] = new Color(0.0f, 0.0f, 1.0f);
                    gl_globals::colors[3] = new Color(1.0f, 1.0f, 0.0f);
                    gl_globals::colors[4] = new Color(1.0f, 0.0f, 1.0f);
                    gl_globals::colors[5] = new Color(0.0f, 1.0f, 1.0f);
                    gl_globals::colors[6] = new Color(1.0f, 1.0f, 1.0f);
                 glClearColor(1.0, 1.0, 1.0, 0.0);
                glShadeModel(GL_SMOOTH);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f,(GLfloat)gl_globals::screen_width/(GLfloat)gl_globals::screen_height,1.0f,1000.0f);
                gluLookAt(0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0);
                glEnable(GL_DEPTH_TEST);
                glEnable(GL_LIGHTING);
                glEnable(GL_POLYGON_SMOOTH);
                glEnable(GL_LINE_SMOOTH);
                glEnable(GL_LIGHT0);
                float pos[] = {100.0, -100.0, -100.0, 1};
                glLightfv(GL_LIGHT0, GL_POSITION, pos); 
                
                glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
            } 

            static void resize (int width, int height)
             {
                gl_globals::screen_width=width; 
                gl_globals::screen_height=height; 
                glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glShadeModel(GL_SMOOTH);
                glViewport(0,0,gl_globals::screen_width,gl_globals::screen_height);
                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                gluPerspective(45.0f,(GLfloat)gl_globals::screen_width/(GLfloat)gl_globals::screen_height,1.0f,1000.0f);
                gluLookAt(0.0, 0.0, -32.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0);
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
                            gl_globals::translation_y_increment = gl_globals::translation_y_increment -0.1;
                        break;

                    case '.':
                            gl_globals::translation_y_increment = gl_globals::translation_y_increment +0.1;
                        break;

                    case 'v':
                            gl_globals::translation_z_increment = gl_globals::translation_z_increment -0.1;
                        break;

                    case 'b':
                            gl_globals::translation_z_increment = gl_globals::translation_z_increment +0.1;
                        break;
                    case 'c':
                            gl_globals::calculate = true;
                        break;
                    case 'p':
                            gl_globals::calculate = false;
                        break;
                    case 's':
                            getAnimator()->rewind();
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
            
            static inline EvolvingAnimator<Tag_, DataType_> * getAnimator()
            {
                return reinterpret_cast<EvolvingAnimator<Tag_, DataType_> *> (gl_globals::animator);
            }

            static void display(void)
            {
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();

                // glTranslatef(-2.0,-2.0,-20);
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
                glRotatef(-45.0f, 0.0, 1.0, 0.0);
                //glRotatef(45.0f, 0.0, 0.0, 1.0);
                gl_globals::translation_x = gl_globals::translation_x + gl_globals::translation_x_increment;
                gl_globals::translation_y = gl_globals::translation_y + gl_globals::translation_y_increment;
                gl_globals::translation_z = gl_globals::translation_z + gl_globals::translation_z_increment;

                glTranslatef(gl_globals::translation_x, 0.0, 0.0);
                glTranslatef(0.0, gl_globals::translation_y, 0.0);
                glTranslatef(0.0, 0.0 , gl_globals::translation_z);
                
                
                //int timeslice = (int)gl_globals::time;     

                glBegin(GL_LINES);
                
                
                EvolvingAnimator<Tag_, DataType_> * animator(getAnimator());
                
                
                for (typename MutableMatrix<DataType_>::ElementIterator i(animator->edges()->begin_non_zero_elements()), i_end(animator->edges()->end_non_zero_elements()); i != i_end ; ++i)
                {
                    if (i.row() > i.column())
                    {
                       glLineWidth((GLfloat)*i * 8);
                       glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, gl_globals::edgeMaterial);
                            
                        DenseVectorRange<DataType_> v1((animator->coordinates())[i.row()]);
                        DenseVectorRange<DataType_> v2((animator->coordinates())[i.column()]);
                        glVertex3f((GLfloat) v1[0],(GLfloat)v1[1], gl_globals::ebene_z);
                        glVertex3f((GLfloat) v2[0],(GLfloat)v2[1], gl_globals::ebene_z);
                    }
                     /*   glLineWidth((GLfloat)*i * 5);
                        glColor3f(1.0, 1.0, 1.0);
                        DenseVectorRange<DataType_> v1(animator->coordinates()[i.row()]);
                        DenseVectorRange<DataType_> v2(animator->coordinates()[i.column()]);
                        glVertex3f((GLfloat) v1[0], (GLfloat)v1[1], gl_globals::ebene_z);
                        glVertex3f((GLfloat) v2[0], (GLfloat)v2[1], gl_globals::ebene_z);*/
                }
                
                glEnd();
                
                // build geometry
                for(unsigned int i = 0; i < animator->coordinates().rows(); ++i)
                {
                    DenseVectorRange<DataType_> dv((animator->coordinates())[i]);
                    glPushMatrix();
                    Color * c =  gl_globals::colors[0];
                    float f[4];
                    f[0] = c->r;
                    f[1] = c->g;
                    f[2] = c->b;
                    f[3] = 1.0f;
                    float shiny[4] = {0.3, 0.3, 0.3, 1.0 };
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, f);
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shiny);
                    //glColor4f();
                    glTranslatef((GLfloat)dv[0], (GLfloat)dv[1], gl_globals::ebene_z);
                    GLUquadricObj  * quad = gluNewQuadric();
                    gluSphere(quad, (GLfloat)(*animator->node_weights())[i] / 10, 8, 8);
                    glPopMatrix();
                }
                
                
                
                
                glFlush();
                glutSwapBuffers();
                std::cout << "time = " << animator->time() << "\n";
                if (gl_globals::calculate)
                {
                    animator->next_step();
                }
                if (animator->end())
                {
                    gl_globals::calculate = false;
                }
            }
    };
}
#endif
