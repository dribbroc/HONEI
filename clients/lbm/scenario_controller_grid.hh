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

#ifndef LBM_SCENARIO_CONTROLLER_GRID_HH
#define LBM_SCENARIO_CONTROLLER_GRID_HH

#include <GL/glut.h>
#include <honei/swe/volume.hh>
#include <honei/lbm/solver_labswe_grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/partial_derivative.hh>
#include <clients/lbm/scenario_controller_base.hh>

template<typename Tag_, typename Prec_> class ScenarioControllerGrid :
    public ScenarioControllerBase
{
    private:
        int scenario_id;


        Grid<D2Q9, Prec_> _grid;
        PackedGridData<D2Q9, Prec_> _data;
        PackedGridInfo<D2Q9> _info;

        unsigned long _dwidth, _dheight, _timestep;
        DenseMatrix<Prec_>* _h;
        DenseMatrix<Prec_>* _u;
        DenseMatrix<Prec_>* _v;
        DenseMatrix<Prec_>* _b;
        DenseMatrix<Prec_>* _b_x;
        DenseMatrix<Prec_>* _b_y;

        DenseMatrix<bool>* _obstacles;

        SolverLABSWEGrid<Tag_, Prec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP>* _solver;

        void _update_scenario()
        {
        }

    public:
        ScenarioControllerGrid(int scen_id) :
            scenario_id(scen_id),
                _b(0),
                _b_x(0),
                _b_y(0)
    {
        srand(time(NULL));
    }
        ~ScenarioControllerGrid()
        {
            delete _h;
            delete _b;
            delete _b_x;
            delete _b_y;
            delete _u;
            delete _v;
            delete _obstacles;
            delete _solver;
        }

        static int get_precision(int scen_id)
        {
            return 0; // todo return the correct accuracy (0(float) or 1(double))
        }

        void init(void)
        {
            delete _b;
            delete _b_x;
            delete _b_y;
            _b = 0;
            _b_x = 0;
            _b_y = 0;

            _grid.destroy();
            _data.destroy();
            _info.destroy();
            _timestep = 0;
            switch (scenario_id)
            {
                case 100:
                    {
                        glutSetWindowTitle("Grid: Laminar flow: Circular dam break 50x50");
                        _dheight = 50;
                        _dwidth = 50;
                        _h = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
                        Cylinder<Prec_> c1(*_h, Prec_(0.06), 25, 25);
                        c1.value();

                        _u = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _v = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b_x = new DenseMatrix<Prec_>(PartialDerivative<Tag_, X, CENTRALDIFF>::value(*_b , Prec_(1)));
                        _b_y = new DenseMatrix<Prec_>(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(*_b , Prec_(1)));

                        _obstacles = new DenseMatrix<bool>(_dheight, _dwidth, false);
                        /*Cylinder<bool> c2(*_obstacles, 1, 6, 10);
                          c2.value();*/

                        /*Cuboid<bool> q2(*_obstacles, 15, 5, 1, 10, 0);
                        q2.value();
                        Cuboid<bool> q3(*_obstacles, 40, 5, 1, 10, 30);
                        q3.value();
                        */
                        _grid.obstacles = _obstacles;
                        _grid.h = _h;
                        _grid.u = _u;
                        _grid.v = _v;
                        _grid.b_x = _b_x;
                        _grid.b_y = _b_y;

                        GridPacker<D2Q9, NOSLIP, Prec_>::pack(_grid, _info, _data);


                        _solver = new SolverLABSWEGrid<Tag_, Prec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> (&_data, &_info, 1., 1., 1.);

                        _solver->do_preprocessing();


                    }
                    break;


                case 101:
                    {
                        glutSetWindowTitle("Grid: Laminar flow: Circular dam break 50x50 with cuboidal obstacles");
                        _dheight = 50;
                        _dwidth = 50;
                        _h = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
                        Cylinder<Prec_> c1(*_h, Prec_(0.06), 25, 25);
                        c1.value();

                        _u = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _v = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b_x = new DenseMatrix<Prec_>(PartialDerivative<Tag_, X, CENTRALDIFF>::value(*_b , Prec_(1)));
                        _b_y = new DenseMatrix<Prec_>(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(*_b , Prec_(1)));

                        _obstacles = new DenseMatrix<bool>(_dheight, _dwidth, false);
                        /*Cylinder<bool> c2(*_obstacles, 1, 6, 10);
                          c2.value();*/

                        Cuboid<bool> q2(*_obstacles, 15, 5, 1, 10, 0);
                        q2.value();
                        Cuboid<bool> q3(*_obstacles, 40, 5, 1, 10, 30);
                        q3.value();
                        _grid.obstacles = _obstacles;
                        _grid.h = _h;
                        _grid.u = _u;
                        _grid.v = _v;
                        _grid.b_x = _b_x;
                        _grid.b_y = _b_y;

                        GridPacker<D2Q9, NOSLIP, Prec_>::pack(_grid, _info, _data);


                        _solver = new SolverLABSWEGrid<Tag_, Prec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> (&_data, &_info, 1., 1., 1.);

                        _solver->do_preprocessing();


                    }
                    break;

                case 102:
                    {
                        glutSetWindowTitle("Grid: Laminar flow: Partial dam break 50x50");
                        _dheight = 50;
                        _dwidth = 50;
                        _h = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
                        /*Cylinder<Prec_> c1(*_h, Prec_(0.06), 25, 25);
                        c1.value();
                        */

                        Cuboid<Prec_> reservoir(*_h, 50, 18, Prec_(0.035), 32, 0);
                        reservoir.value();

                        _u = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _v = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));

                        _b_x = new DenseMatrix<Prec_>(PartialDerivative<Tag_, X, CENTRALDIFF>::value(*_b , Prec_(1)));
                        _b_y = new DenseMatrix<Prec_>(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(*_b , Prec_(1)));
                        _obstacles = new DenseMatrix<bool>(_dheight, _dwidth, false);
                        /*Cylinder<bool> c2(*_obstacles, 1, 6, 10);
                          c2.value();*/

                        Cuboid<bool> q2(*_obstacles, 15, 2, 1, 30, 0);
                        q2.value();
                        Cuboid<bool> q3(*_obstacles, 40, 2, 1, 30, 20);
                        q3.value();
                        _grid.obstacles = _obstacles;
                        _grid.h = _h;
                        _grid.u = _u;
                        _grid.v = _v;
                        _grid.b_x = _b_x;
                        _grid.b_y = _b_y;

                        GridPacker<D2Q9, NOSLIP, Prec_>::pack(_grid, _info, _data);


                        _solver = new SolverLABSWEGrid<Tag_, Prec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> (&_data, &_info, 1., 1., 1.);

                        _solver->do_preprocessing();


                    }
                    break;

                case 103:
                    {
                        glutSetWindowTitle("Grid: Laminar flow: Circular dam break over uneven bed 50x50");
                        _dheight = 50;
                        _dwidth = 50;
                        _h = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
                        Cylinder<Prec_> c1(*_h, Prec_(0.06), 25, 25);
                        c1.value();

                        _u = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _v = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        _b = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.));
                        srand(time(NULL));
                        for (unsigned long i(0); i < _dheight; ++i)
                        {
                            for (unsigned long j(0); j < _dwidth; ++j)
                            {
                                (*_b)( i , j) += (Prec_)(Prec_(rand()) / RAND_MAX * 0.003);
                            }
                        }
                        _b_x = new DenseMatrix<Prec_>(PartialDerivative<Tag_, X, CENTRALDIFF>::value(*_b , Prec_(1)));
                        _b_y = new DenseMatrix<Prec_>(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(*_b , Prec_(1)));

                        _obstacles = new DenseMatrix<bool>(_dheight, _dwidth, false);
                        /*Cylinder<bool> c2(*_obstacles, 1, 6, 10);
                          c2.value();*/

                        /*Cuboid<bool> q2(*_obstacles, 15, 5, 1, 10, 0);
                        q2.value();
                        Cuboid<bool> q3(*_obstacles, 40, 5, 1, 10, 30);
                        q3.value();
                        */
                        _grid.obstacles = _obstacles;
                        _grid.h = _h;
                        _grid.u = _u;
                        _grid.v = _v;
                        _grid.b_x = _b_x;
                        _grid.b_y = _b_y;

                        GridPacker<D2Q9, NOSLIP, Prec_>::pack(_grid, _info, _data);


                        _solver = new SolverLABSWEGrid<Tag_, Prec_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> (&_data, &_info, 1., 1., 1.);

                        _solver->do_preprocessing();


                    }
                    break;
            }

        }


        void do_timestep(void)
        {
            _update_scenario();
            _solver->solve();
            GridPacker<D2Q9, NOSLIP, Prec_>::unpack(_grid, _info, _data);
            ++_timestep;
        }

        void render(bool show_ground, bool use_quads, bool enable_alpha_blending, bool show_water, float alpha)
        {

            glScalef(1.0f, 1.0f, 100.0f);
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
                            glVertex3d(i,j,(*_b)[j][i]);
                            glColor3f(1.0, 0.8, 0.0);
                            glVertex3d(i+1,j, (*_b)[j][i+1]);
                            glVertex3d(i+1,j+1, (*_b)[j+1][i+1]);
                            glVertex3d(i,j+1, (*_b)[j+1][i]);
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
                            glVertex3d(i,j, (*_b)[j][i]);
                            glColor3f(1.0, 0.0, 0.0);
                            glVertex3d(i+1,j, (*_b)[j][i+1]);
                        }
                        ++i;
                        if (i >= _dwidth-1)
                            break;
                        for(int j2 = _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, (*_b)[j2][i]);
                            glColor3f(1.0, 0.8, 0.0);
                            glVertex3d(i+1,j2, (*_b)[j2][i+1]);
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
                            glVertex3d(i,j,((*_h)(j,i) + (*_b)(j, i)));
                            glColor4f(0.0, 1.0, 1.0, alpha);
                            glVertex3d(i+1,j,((*_h)[j][i+1] + (*_b)[j][i+1]));
                            glVertex3d(i+1,j+1,((*_h)[j+1][i+1] + (*_b)[j+1][i+1]));
                            glVertex3d(i,j+1,((*_h)[j+1][i] + (*_b)[j+1][i]));
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
                            glVertex3d(i,j, ((*_h)[j][i] +  (*_b)[j][i]));
                            glColor4f(0.0, 0.0, 1.0,  alpha);
                            glVertex3d(i+1,j,((*_h)[j][i+1] +  (*_b)[j][i+1]));
                        }
                        ++i;
                        if (i >=  _dwidth-1)
                            break;
                        for(int j2 =  _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, ((*_h)[j2][i] + (*_b)[j2][i]));
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i+1,j2, ((*_h)[j2][i+1] +  (*_b)[j2][i+1]));
                        }
                    }
                    glEnd();
                }
            }
            glBegin(GL_QUADS);
            for(unsigned int i = 0; i < _dwidth-1; ++i)
            {
                for(unsigned int j = 0; j < _dheight-1; ++j)
                {
                    glColor3f(1.0, 0.0, 0.0);
                    if ((*_obstacles)(j, i))
                        glVertex3d(i,j, 0.1);
                    else
                        glVertex3d(i,j,0);
                    glColor3f(1.0, 0.8, 0.0);
                    if ((*_obstacles)(j, i + 1))
                        glVertex3d(i+1,j, 0.1);
                    else
                        glVertex3d(i+1,j, 0);
                    if ((*_obstacles)(j + 1, i + 1))
                        glVertex3d(i+1,j+1, 0.1);
                    else
                        glVertex3d(i+1,j+1, 0);
                    if ((*_obstacles)(j + 1, i))
                        glVertex3d(i,j+1, 0.1);
                    else
                        glVertex3d(i,j+1, 0);
                }
            }
            glEnd();
        }
};
#endif
