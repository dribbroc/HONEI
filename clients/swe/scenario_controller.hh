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

        unsigned long _dwidth, _dheight;

        Prec_ _dt, _dx, _dy, _manning;

        double _eps;

        RelaxSolver<Tag_, Prec_, Prec_, Prec_, Prec_, Prec_, source_types::SIMPLE, boundaries::REFLECT, FIXED>* _solver;

    public:
        ScenarioController(int scen_id) :
            scenario_id(scen_id)
    {
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
            return scen_id;
        }

        void init(void)
        {
            switch (scenario_id)
            {
                //Rain 90x90:
                case 0:
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
                    unsigned long entries = 3*((_dwidth*_dheight)+4*(_dwidth+_dheight+4));
                    _eps = 10e-6;
                    _manning = 0;

                    _u = new DenseVector<Prec_>(entries, Prec_(0));
                    _v = new DenseVector<Prec_>(entries, Prec_(0));
                    _w = new DenseVector<Prec_>(entries, Prec_(0));

                    _bx = new DenseVector<Prec_>(entries/3, Prec_(0));
                    _by = new DenseVector<Prec_>(entries/3, Prec_(0));

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
        }

        void do_timestep(void)
        {
            _solver->solve();
        }

        DenseMatrix<Prec_>* get_field_water()
        {
            return _height;
        }

        DenseMatrix<Prec_>* get_field_bottom()
        {
            return _bottom;
        }
};
#endif
