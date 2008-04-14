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

#ifndef SWE_SCENARIO_CONTROLLER_HH
#define SWE_SCENARIO_CONTROLLER_HH

#include <GL/glut.h>
#include <honei/libswe/volume.hh>
#include <honei/libswe/relax_solver.hh>
#include <honei/libswe/scenario_manager.hh>

template<typename Tag_, typename Prec_> class ScenarioController
{
    private:
        int scenario_id;

        Scenario<Prec_, RELAX, REFLECT> * _scenario;

        DenseMatrix<Prec_>* _height;
        DenseMatrix<Prec_>* _bottom;
        DenseMatrix<Prec_>* _u1;
        DenseMatrix<Prec_>* _u2;

        DenseVector<Prec_>* _u;
        DenseVector<Prec_>* _v;
        DenseVector<Prec_>* _w;
        DenseVector<Prec_>* _bx;
        DenseVector<Prec_>* _by;
        DenseVector<Prec_>* _c;
        DenseVector<Prec_>* _d;

        unsigned long _dwidth, _dheight, _entries, _timestep;

        Prec_ _dt, _dx, _dy, _manning;

        double _eps;

        RelaxSolver<Tag_, Prec_, Prec_, Prec_, Prec_, Prec_, source_types::SIMPLE, boundaries::REFLECT, FIXED>* _solver;

        void _update_scenario()
        {
            switch(scenario_id)
            {
                case 1:
                    {
                        if(_timestep % 25==0)
                        {
                            //generate numbers in [6,_d_width-6]
                            Prec_ x = rand() % (_dwidth - 14);
                            Prec_ y = rand() % (_dheight - 14);

                            Prec_ strength = rand() % 20;
                            Cylinder<Prec_> c(*_height, (unsigned long)strength, (unsigned long)(x + 7), (unsigned long)(y +7));
                            c.value();
                            _solver->do_preprocessing();

                        }
                    }
                    break;
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
            delete _u1;
            delete _u2;
            delete _u;
            delete _v;
            delete _w;
            delete _bx;
            delete _by;
            delete _c;
            delete _d;
            delete _solver;
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
                        _dwidth = 90;
                        _dheight = 90;
                        _dt = 5./24.;
                        _dx = 5;
                        _dy = 5;

                        _c = new DenseVector<Prec_>(3);
                        (*_c)[0] = 12.;
                        (*_c)[1] = 7.;
                        (*_c)[2] = 12.;
                        _d = new DenseVector<Prec_>(_c->copy());

                        _height = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                        _bottom = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                        _u1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _u2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _entries = 3*((_dwidth*_dheight)+4*(_dwidth+_dheight+4));
                        _eps = 10e-6;
                        _manning = 0;

                        _u = new DenseVector<Prec_>(_entries, Prec_(0));
                        _v = new DenseVector<Prec_>(_entries, Prec_(0));
                        _w = new DenseVector<Prec_>(_entries, Prec_(0));

                        _bx = new DenseVector<Prec_>(_entries/3, Prec_(0));
                        _by = new DenseVector<Prec_>(_entries/3, Prec_(0));

                        Cylinder<Prec_> c1(*_height, Prec_(15.), _dwidth/2, _dheight/2);
                        c1.value();

                        ScenarioManager<Prec_, swe_solvers::RELAX, boundaries::REFLECT> scen_man;
                        _scenario =  new Scenario<Prec_, RELAX, REFLECT>(_dwidth, _dheight);

                        scen_man.allocate_scenario(_scenario);
                        scen_man.allocate_scalarfields(_height, _bottom, _u1, _u2);
                        scen_man.allocate_relax_vectors(_u, _v, _w, _c, _d);
                        scen_man.allocate_bottom_slopes(_bx, _by);
                        scen_man.set_environmental_variables(_dx, _dy, _dt, _manning, _eps);

                        _solver = new RelaxSolver<Tag_, Prec_, Prec_, Prec_, Prec_, Prec_, source_types::SIMPLE, boundaries::REFLECT, FIXED>(*_scenario);

                        if(scen_man.validate())
                        {
                            _solver->do_preprocessing();
                        }
                    }
                    break;
                    //Rain 64x64:
                case 1:
                    {
                        _dwidth = 64;
                        _dheight = 64;
                        _dt = 5./24.;
                        _dx = 5;
                        _dy = 5;

                        _c = new DenseVector<Prec_>(3);
                    (*_c)[0] = 12.;
                    (*_c)[1] = 7.;
                    (*_c)[2] = 12.;
                    _d = new DenseVector<Prec_>(_c->copy());

                    _height = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                    _bottom = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                    _u1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                    _u2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                    _entries = 3*((_dwidth*_dheight)+4*(_dwidth+_dheight+4));
                    _eps = 10e-6;
                    _manning = 0;

                    _u = new DenseVector<Prec_>(_entries, Prec_(0));
                    _v = new DenseVector<Prec_>(_entries, Prec_(0));
                    _w = new DenseVector<Prec_>(_entries, Prec_(0));

                    _bx = new DenseVector<Prec_>(_entries/3, Prec_(0));
                    _by = new DenseVector<Prec_>(_entries/3, Prec_(0));

                    Cylinder<Prec_> c1(*_height, Prec_(15.), _dwidth/2, _dheight/2);
                    c1.value();

                    ScenarioManager<Prec_, swe_solvers::RELAX, boundaries::REFLECT> scen_man;
                    _scenario =  new Scenario<Prec_, RELAX, REFLECT>(_dwidth, _dheight);

                    scen_man.allocate_scenario(_scenario);
                    scen_man.allocate_scalarfields(_height, _bottom, _u1, _u2);
                    scen_man.allocate_relax_vectors(_u, _v, _w, _c, _d);
                    scen_man.allocate_bottom_slopes(_bx, _by);
                    scen_man.set_environmental_variables(_dx, _dy, _dt, _manning, _eps);

                    _solver = new RelaxSolver<Tag_, Prec_, Prec_, Prec_, Prec_, Prec_, source_types::SIMPLE, boundaries::REFLECT, FIXED>(*_scenario);

                    if(scen_man.validate())
                    {
                        _solver->do_preprocessing();
                    }
                    }
                    break;

                //Full dam break 90x90:
                case 2:
                    {
                        _dwidth = 90;
                        _dheight = 90;
                        _dt = 5./24.;
                        _dx = 5;
                        _dy = 5;

                        _c = new DenseVector<Prec_>(3);
                        (*_c)[0] = 12.;
                        (*_c)[1] = 7.;
                        (*_c)[2] = 12.;
                        _d = new DenseVector<Prec_>(_c->copy());

                        _height = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                        _bottom = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(5.));
                        _u1 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _u2 = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _entries = 3*((_dwidth*_dheight)+4*(_dwidth+_dheight+4));
                        _eps = 10e-6;
                        _manning = 0;

                        _u = new DenseVector<Prec_>(_entries, Prec_(0));
                        _v = new DenseVector<Prec_>(_entries, Prec_(0));
                        _w = new DenseVector<Prec_>(_entries, Prec_(0));

                        _bx = new DenseVector<Prec_>(_entries/3, Prec_(0));
                        _by = new DenseVector<Prec_>(_entries/3, Prec_(0));

                        Cuboid<Prec_> q2(*_height, Prec_(15), Prec_(88), Prec_(15),1,89);
                        q2.value();

                        ScenarioManager<Prec_, swe_solvers::RELAX, boundaries::REFLECT> scen_man;
                        _scenario =  new Scenario<Prec_, RELAX, REFLECT>(_dwidth, _dheight);

                        scen_man.allocate_scenario(_scenario);
                        scen_man.allocate_scalarfields(_height, _bottom, _u1, _u2);
                        scen_man.allocate_relax_vectors(_u, _v, _w, _c, _d);
                        scen_man.allocate_bottom_slopes(_bx, _by);
                        scen_man.set_environmental_variables(_dx, _dy, _dt, _manning, _eps);

                        _solver = new RelaxSolver<Tag_, Prec_, Prec_, Prec_, Prec_, Prec_, source_types::SIMPLE, boundaries::REFLECT, FIXED>(*_scenario);

                        if(scen_man.validate())
                        {
                            _solver->do_preprocessing();
                        }
                    }
                    break;


            }

        }


void do_timestep(void)
{
    _update_scenario();
    _solver->solve();
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
