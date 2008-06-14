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

        DenseMatrix<Prec_>* _d_0;
        DenseMatrix<Prec_>* _d_1;
        DenseMatrix<Prec_>* _d_2;
        DenseMatrix<Prec_>* _d_3;
        DenseMatrix<Prec_>* _d_4;
        DenseMatrix<Prec_>* _d_5;
        DenseMatrix<Prec_>* _d_6;
        DenseMatrix<Prec_>* _d_7;
        DenseMatrix<Prec_>* _d_8;

        DenseMatrix<Prec_>* _e_d_0;
        DenseMatrix<Prec_>* _e_d_1;
        DenseMatrix<Prec_>* _e_d_2;
        DenseMatrix<Prec_>* _e_d_3;
        DenseMatrix<Prec_>* _e_d_4;
        DenseMatrix<Prec_>* _e_d_5;
        DenseMatrix<Prec_>* _e_d_6;
        DenseMatrix<Prec_>* _e_d_7;
        DenseMatrix<Prec_>* _e_d_8;

        DenseMatrix<Prec_>* _t_d_0;
        DenseMatrix<Prec_>* _t_d_1;
        DenseMatrix<Prec_>* _t_d_2;
        DenseMatrix<Prec_>* _t_d_3;
        DenseMatrix<Prec_>* _t_d_4;
        DenseMatrix<Prec_>* _t_d_5;
        DenseMatrix<Prec_>* _t_d_6;
        DenseMatrix<Prec_>* _t_d_7;
        DenseMatrix<Prec_>* _t_d_8;

        DenseVector<Prec_>* _v_x;
        DenseVector<Prec_>* _v_y;

        DenseMatrix<Prec_>* _s_x;
        DenseMatrix<Prec_>* _s_y;
        DenseMatrix<Prec_>* _d_x;
        DenseMatrix<Prec_>* _d_y;

        unsigned long _dwidth, _dheight, _timestep;

        SolverLABSWE<Tag_, Prec_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> * _solver_1;

        void _update_scenario()
        {
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
            delete _d_0;
            delete _d_1;
            delete _d_2;
            delete _d_3;
            delete _d_4;
            delete _d_5;
            delete _d_6;
            delete _d_7;
            delete _d_8;
            delete _e_d_0;
            delete _e_d_1;
            delete _e_d_2;
            delete _e_d_3;
            delete _e_d_4;
            delete _e_d_5;
            delete _e_d_6;
            delete _e_d_7;
            delete _e_d_8;
            delete _t_d_0;
            delete _t_d_1;
            delete _t_d_2;
            delete _t_d_3;
            delete _t_d_4;
            delete _t_d_5;
            delete _t_d_6;
            delete _t_d_7;
            delete _t_d_8;
            delete _v_x;
            delete _v_y;
            delete _s_x;
            delete _s_y;
            delete _d_x;
            delete _d_y;
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
                case 0:
                    {
                        _dheight = 50;
                        _dwidth = 50;
                        _height = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
                        Cylinder<Prec_> c1(*_height, Prec_(0.075), 25, 25);
                        c1.value();

                        _bottom = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _u = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _v = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        //All needed distribution functions:

                        _d_0 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_3 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_4 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_5 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_6 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_7 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_8 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        _e_d_0 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_3 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_4 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_5 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_6 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_7 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _e_d_8 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        _t_d_0 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_3 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_4 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_5 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_6 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_7 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _t_d_8 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        //All needed vectors:
                        _v_x = new DenseVector<Prec_>(9, Prec_(0));
                        _v_y = new DenseVector<Prec_>(9, Prec_(0));

                        //Other matrices needed by solver:

                        _s_x = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _s_y = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_x = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _d_y = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        _solver_1 = new SolverLABSWE<Tag_, Prec_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>(1.,1.,1., _dwidth, _dheight, _height, _bottom, _u, _v);
                        _solver_1->set_distribution(_d_0, _d_1, _d_2, _d_3, _d_4, _d_5, _d_6, _d_7, _d_8);
                        _solver_1->set_eq_distribution(_e_d_0, _e_d_1, _e_d_2, _e_d_3, _e_d_4, _e_d_5, _e_d_6, _e_d_7, _e_d_8);
                        _solver_1->set_temp_distribution(_t_d_0, _t_d_1, _t_d_2, _t_d_3, _t_d_4, _t_d_5, _t_d_6, _t_d_7, _t_d_8);
                        _solver_1->set_vectors(_v_x, _v_y);
                        _solver_1->set_source(_s_x, _s_y);
                        _solver_1->set_slopes(_d_x, _d_y);
                        _solver_1->do_preprocessing();
                    }
                    break;
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
                            glVertex3d(i,j, 100.*((*_height)[j][i] + (*_bottom)[j][i]));
                            glColor4f(0.0, 1.0, 1.0, alpha);
                            glVertex3d(i+1,j, 100.*((*_height)[j][i+1] + (*_bottom)[j][i+1]));
                            glVertex3d(i+1,j+1, 100.*((*_height)[j+1][i+1] + (*_bottom)[j+1][i+1]));
                            glVertex3d(i,j+1, 100.*((*_height)[j+1][i] + (*_bottom)[j+1][i]));
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
                            glVertex3d(i,j, 100.*((*_height)[j][i] +  (*_bottom)[j][i]));
                            glColor4f(0.0, 0.0, 1.0,  alpha);
                            glVertex3d(i+1,j, 100.*((*_height)[j][i+1] +  (*_bottom)[j][i+1]));
                        }
                        ++i;
                        if (i >=  _dwidth-1)
                            break;
                        for(int j2 =  _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, 100.*((*_height)[j2][i] +  (*_bottom)[j2][i]));
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i+1,j2, 100.*((*_height)[j2][i+1] +  (*_bottom)[j2][i+1]));
                        }
                    }
                    glEnd();
                }
            }
        }
};
#endif
