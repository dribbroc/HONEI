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

#ifndef LBM_SCENARIO_CONTROLLER_HH
#define LBM_SCENARIO_CONTROLLER_HH

#include <GL/glut.h>
#include <honei/swe/volume.hh>
#include <honei/lbm/solver_labswe.hh>

template<typename Tag_, typename Prec_> class ScenarioController
{
    private:
        int scenario_id;

        DenseMatrix<Prec_>* _height;
        DenseMatrix<Prec_>* _bottom;
        DenseMatrix<Prec_>* _u;
        DenseMatrix<Prec_>* _v;

        unsigned long _dwidth, _dheight, _timestep;

        SolverLABSWE<Tag_, Prec_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> * _solver_1;

            void _update_scenario()
            {
                switch(scenario_id)
                {
                    case 1:
                        {
                        }
                }
            }

    public:
        ScenarioController(int scen_id) :
            scenario_id(scen_id)
    {
        srand(time(NULL));
    }
        ~ScenarioController()
        {
            delete _height;
            delete _bottom;
            delete _u;
            delete _v;
            delete _solver_1;
        }

        static int get_precision(int scen_id)
        {
            return 0; // todo return the correct accuracy (0(float) or 1(double))
        }

        void init(void)
        {
            _timestep = 0;
            //todo delete old data
            switch (scenario_id)
            {
                //Rain 90x90:
                case 0:
                    {
                        //TODO: Preprocessing, setup, etc
                        _solver_1 = new SolverLABSWE<Tag_, Prec_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC>(1.,1.,1., _dwidth, _dheight, _height, _bottom, _u, _v);
                    }
            }

        }


        void do_timestep(void)
        {
            _update_scenario();
            _solver_1->solve();
            ++_timestep;
        }

        void render(bool show_ground, bool use_quads, bool enable_alpha_blending, bool show_water, float alpha)
        {
            if (show_ground)
            {
                if(use_quads)
                {
                    glBegin(GL_QUADS);
                    for(unsigned int i = 0; i < _dwidth-1; ++i)
                    {
                        for(unsigned int j = 0; j < _dheight-1; ++j)
                        {
                            glColor3f(1.0, 0.0, 0.0);
                            glVertex3d(i,j,(*_bottom)[j][i]);
                            glColor3f(1.0, 0.8, 0.0);
                            glVertex3d(i+1,j, (*_bottom)[j][i+1]);
                            glVertex3d(i+1,j+1, (*_bottom)[j+1][i+1]);
                            glVertex3d(i,j+1, (*_bottom)[j+1][i]);
                        }
                    }
                    glEnd();
                }
                else
                {
                    glBegin(GL_TRIANGLE_STRIP);
                    for(unsigned int i = 0; i < _dwidth-1; ++i)
                    {
                        for(unsigned int j = 0; j < _dheight; j++)
                        {
                            glColor3f(1.0, 0.8, 0.0);
                            glVertex3d(i,j, (*_bottom)[j][i]);
                            glColor3f(1.0, 0.0, 0.0);
                            glVertex3d(i+1,j, (*_bottom)[j][i+1]);
                        }
                        ++i;
                        if (i >= _dwidth-1)
                            break;
                        for(int j2 = _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, (*_bottom)[j2][i]);
                            glColor3f(1.0, 0.8, 0.0);
                            glVertex3d(i+1,j2, (*_bottom)[j2][i+1]);
                        }
                    }
                    glEnd();
                }
            }
            if(enable_alpha_blending)
            {
                glEnable (GL_BLEND);
                glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            }
            else
                glDisable (GL_BLEND);

            if (show_water)
            {
                if(use_quads)
                {
                    glBegin(GL_QUADS);
                    for(unsigned int i = 0; i < _dwidth-1; ++i)
                    {
                        for(unsigned int j = 0; j <_dheight-1; ++j)
                        {
                            glColor4f(0.0, 0.0, 1.0, alpha);
                            glVertex3d(i,j, (*_height)[j][i] + (*_bottom)[j][i]);
                            glColor4f(0.0, 1.0, 1.0, alpha);
                            glVertex3d(i+1,j, (*_height)[j][i+1] + (*_bottom)[j][i+1]);
                            glVertex3d(i+1,j+1, (*_height)[j+1][i+1] + (*_bottom)[j+1][i+1]);
                            glVertex3d(i,j+1, (*_height)[j+1][i] + (*_bottom)[j+1][i]);
                        }
                    }
                    glEnd();
                }
                else
                {
                    glBegin(GL_TRIANGLE_STRIP);
                    for(unsigned int i = 0; i <  _dwidth-1; ++i)
                    {
                        for(unsigned int j = 0; j <  _dheight; j++)
                        {
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i,j, (*_height)[j][i] +  (*_bottom)[j][i]);
                            glColor4f(0.0, 0.0, 1.0,  alpha);
                            glVertex3d(i+1,j, (*_height)[j][i+1] +  (*_bottom)[j][i+1]);
                        }
                        ++i;
                        if (i >=  _dwidth-1)
                            break;
                        for(int j2 =  _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, (*_height)[j2][i] +  (*_bottom)[j2][i]);
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i+1,j2, (*_height)[j2][i+1] +  (*_bottom)[j2][i+1]);
                        }
                    }
                    glEnd();
                }
            }
        }

};
#endif
