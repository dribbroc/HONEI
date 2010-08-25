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
#include <honei/math/multigrid.hh>
#include <honei/math/fill_matrix.hh>
#include <honei/math/fill_vector.hh>
#include <clients/poisson/scenario_controller_base.hh>
#include <iostream>
template<typename Tag_, typename Prec_> class ScenarioController :
    public ScenarioControllerBase
{
    private:
        int scenario_id;

        DenseMatrix<float> * _u;
        unsigned long _timestep;

        unsigned long _root_n;

        void _update_scenario()
        {
        }

    public:
        ScenarioController(int scen_id) :
            scenario_id(scen_id)
    {
        srand(time(NULL));
    }
        virtual ~ScenarioController()
        {
            delete _u;
        }

        static int get_precision(int /*scen_id*/)
        {
            return 1; // todo return the correct accuracy (0(float) or 1(double))
        }

        void init(void)
        {
            _timestep = 0;
            //todo delete old data
            switch (scenario_id)
            {
                case 100:
                    {
                        glutSetWindowTitle("F(x,y) = f, Dirichlet 2, Neumann east MG V2+2, N=33^2");
                        _root_n = 33;
                        unsigned long n(_root_n * _root_n);
                        _u = new DenseMatrix<float>(_root_n, _root_n, Prec_(0.00));
                        MGInfo<Prec_, BandedMatrixQ1<Prec_> > info;
                        //configuration constants: /TODO: set/allocate!!!
                        info.is_smoother = false;
                        DenseVector<unsigned long> mask(8);

                        info.macro_border_mask = new DenseVector<unsigned long>(8);
                        for(unsigned long i(0); i < 8; ++i)
                        {
                            (*info.macro_border_mask)[i] = 2;
                        }
                        //set Neumann boundaries:
                        (*info.macro_border_mask)[5] =1;

                        info.min_level = 1;
                        switch(n)
                        {
                            case 1050625:
                                {
                                    info.max_level = 10;
                                }
                                break;
                            case 263169:
                                {
                                    info.max_level = 9;
                                }
                                break;
                            case 66049:
                                {
                                    info.max_level = 8;
                                }
                                break;
                            case 16641:
                                {
                                    info.max_level = 7;
                                }
                                break;
                            case 4225:
                                {
                                    info.max_level = 6;
                                }
                                break;
                            case 1089:
                                {
                                    info.max_level = 5;
                                }
                                break;
                            case 289:
                                {
                                    info.max_level = 4;
                                }
                                break;
                            case 81:
                                {
                                    info.max_level = 3;
                                }
                                break;
                            case 25:
                                {
                                    info.max_level = 2;
                                }
                                break;
                            case 9:
                                {
                                    info.max_level = 1;
                                }
                                break;
                        }

                        info.n_max_iter = 16;
                        info.initial_zero = false;
                        info.tolerance = 1e-8;
                        info.convergence_check = true;

                        info.n_pre_smooth = 2;
                        info.n_post_smooth = 2;
                        info.n_max_iter_coarse = ((unsigned long)sqrt((Prec_)(pow((Prec_)2 , (Prec_)info.max_level) + 1)*(pow((Prec_)2 , (Prec_)info.max_level) + 1)));
                        info.tolerance_coarse = 1e-2;
                        info.adapt_correction_factor = 1.;

                        for (unsigned long i(0) ; i < info.min_level; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(i == 0)
                                size = 9;

                            DenseVector<Prec_> dummy_band(size, Prec_(0));
                            BandedMatrixQ1<Prec_> ac_a(size, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band);
                            info.a.push_back(ac_a);
                            // iteration vectors
                            DenseVector<Prec_> ac_c(size, Prec_(0));
                            info.c.push_back(ac_c);
                            DenseVector<Prec_> ac_d(size, Prec_(0));
                            info.d.push_back(ac_d);
                            DenseVector<Prec_> ac_rhs(size, Prec_(0));
                            info.rhs.push_back(ac_rhs);
                            DenseVector<Prec_> ac_x(size, Prec_(0));
                            info.x.push_back(ac_x);

                            info.diags_inverted.push_back(dummy_band.copy());
                        }


                        for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
                        {
                            unsigned long size = (unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1));
                            // iteration vectors
                            DenseVector<Prec_> ac_c(size, Prec_(0));
                            info.c.push_back(ac_c);
                            DenseVector<Prec_> ac_d(size, Prec_(0));
                            info.d.push_back(ac_d);
                            DenseVector<Prec_> ac_x(size, Prec_(0));
                            info.x.push_back(ac_x);

                            DenseVector<Prec_> dummy_band(size, Prec_(0));
                            info.diags_inverted.push_back(dummy_band.copy());
                        }

                        //assemble all needed levels' matrices:
                        for(unsigned long i(info.min_level); i <= info.max_level; ++i)
                        {
                            unsigned long N = (unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1));
                            DenseVector<Prec_> band(N);
                            BandedMatrixQ1<Prec_> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                            DenseVector<Prec_> current_rhs(N);


                            FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_matrix);

                            FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_rhs);

                            info.rhs.push_back(current_rhs);
                            info.a.push_back(current_matrix);
                        }
                        //clear x data
                        for(unsigned long i(0) ; i < info.max_level ; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(size==0)
                                size = 9;

                            DenseVector<Prec_> null(size , Prec_(0));
                            info.x[i] = null.copy();
                        }
                        //SET DIAG_INVERTED:
                        for (unsigned long i(0) ; i <= info.max_level; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(i == 0)
                                size = 9;

                            DenseVector<Prec_> scaled_diag_inverted(info.a[i].band(DD).copy());
                            ElementInverse<Tag_>::value(scaled_diag_inverted);
                            Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                            info.diags_inverted[i] = scaled_diag_inverted.copy();
                        }

                        DenseVector<Prec_> result(n, Prec_(0));
                        Multigrid<Tag_, Tag_, NONE, JAC, CYCLE::V, FIXED >::value(info.a[info.max_level], info.rhs[info.max_level], result, (unsigned long)11, std::numeric_limits<Prec_>::epsilon(), info);

                        //write result to scalarfield:
                        unsigned long row_index(0);
                        unsigned long column_index(0);
                        unsigned long index(0);

                        for (; index < info.rhs[info.max_level].size() ; ++index)
                        {
                            (*_u)(row_index, column_index) = result[index] *100;
                            if((index + 1) % _root_n == 0)
                            {
                                ++row_index;
                                column_index = 0;
                            }
                            else
                                ++column_index;
                        }
                    }
                    break;

                case 101:
                    {
                        glutSetWindowTitle("F(x,y) = f, Dirichlet 2, Neumann east MG V2+2, N=33^2, mixedprec");
                        _root_n = 33;
                        unsigned long n(_root_n * _root_n);
                        _u = new DenseMatrix<float>(_root_n, _root_n, Prec_(0.00));
                        MGInfo<float, BandedMatrixQ1<float> > info;
                        //configuration constants: /TODO: set/allocate!!!
                        info.is_smoother = false;
                        DenseVector<unsigned long> mask(8);

                        info.macro_border_mask = new DenseVector<unsigned long>(8);
                        for(unsigned long i(0); i < 8; ++i)
                        {
                            (*info.macro_border_mask)[i] = 2;
                        }
                        //set Neumann boundaries:
                        (*info.macro_border_mask)[5] =1;

                        info.min_level = 1;
                        switch(n)
                        {
                            case 1050625:
                                {
                                    info.max_level = 10;
                                }
                                break;
                            case 263169:
                                {
                                    info.max_level = 9;
                                }
                                break;
                            case 66049:
                                {
                                    info.max_level = 8;
                                }
                                break;
                            case 16641:
                                {
                                    info.max_level = 7;
                                }
                                break;
                            case 4225:
                                {
                                    info.max_level = 6;
                                }
                                break;
                            case 1089:
                                {
                                    info.max_level = 5;
                                }
                                break;
                            case 289:
                                {
                                    info.max_level = 4;
                                }
                                break;
                            case 81:
                                {
                                    info.max_level = 3;
                                }
                                break;
                            case 25:
                                {
                                    info.max_level = 2;
                                }
                                break;
                            case 9:
                                {
                                    info.max_level = 1;
                                }
                                break;
                        }

                        info.n_max_iter = 2;
                        info.initial_zero = false;
                        info.tolerance = 1e-2;
                        info.convergence_check = false;

                        info.n_pre_smooth = 2;
                        info.n_post_smooth = 2;
                        info.n_max_iter_coarse = ((unsigned long)sqrt((Prec_)(pow((Prec_)2 , (Prec_)info.max_level) + 1)*(pow((Prec_)2 , (Prec_)info.max_level) + 1)));
                        info.tolerance_coarse = 1e-2;
                        info.adapt_correction_factor = 1.;

                        for (unsigned long i(0) ; i < info.min_level; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(i == 0)
                                size = 9;

                            DenseVector<float> dummy_band(size, float(0));
                            BandedMatrixQ1<float> ac_a(size, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band, dummy_band);
                            info.a.push_back(ac_a);
                            // iteration vectors
                            DenseVector<float> ac_c(size, float(0));
                            info.c.push_back(ac_c);
                            DenseVector<float> ac_d(size, float(0));
                            info.d.push_back(ac_d);
                            DenseVector<float> ac_rhs(size, float(0));
                            info.rhs.push_back(ac_rhs);
                            DenseVector<float> ac_x(size, float(0));
                            info.x.push_back(ac_x);

                            info.diags_inverted.push_back(dummy_band.copy());
                        }


                        for (unsigned long i(info.min_level) ; i <= info.max_level; ++i)
                        {
                            unsigned long size = (unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1));
                            // iteration vectors
                            DenseVector<float> ac_c(size, float(0));
                            info.c.push_back(ac_c);
                            DenseVector<float> ac_d(size, float(0));
                            info.d.push_back(ac_d);
                            DenseVector<float> ac_x(size, float(0));
                            info.x.push_back(ac_x);

                            DenseVector<float> dummy_band(size, float(0));
                            info.diags_inverted.push_back(dummy_band.copy());
                        }

                        //assemble all needed levels' matrices:
                        for(unsigned long i(info.min_level); i <= info.max_level; ++i)
                        {
                            unsigned long N = (unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1));
                            DenseVector<float> band(N);
                            BandedMatrixQ1<float> current_matrix(N, band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy(), band.copy());
                            DenseVector<float> current_rhs(N);


                            FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_matrix);

                            FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(current_rhs);

                            info.rhs.push_back(current_rhs);
                            info.a.push_back(current_matrix);
                        }
                        //clear x data
                        for(unsigned long i(0) ; i < info.max_level ; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(size==0)
                                size = 9;

                            DenseVector<float> null(size , float(0));
                            info.x[i] = null.copy();
                        }
                        //SET DIAG_INVERTED:
                        for (unsigned long i(0) ; i <= info.max_level; ++i)
                        {
                            unsigned long size((unsigned long)(((unsigned long)pow((Prec_)2, (Prec_)i) + 1) * ((unsigned long)pow((Prec_)2, (Prec_)i) + 1)));
                            if(i == 0)
                                size = 9;

                            DenseVector<float> scaled_diag_inverted(info.a[i].band(DD).copy());
                            ElementInverse<Tag_>::value(scaled_diag_inverted);
                            Scale<Tag_>::value(scaled_diag_inverted, 0.7);

                            info.diags_inverted[i] = scaled_diag_inverted.copy();
                        }

                        DenseVector<Prec_> null(info.rhs[info.max_level].size() , Prec_(0));
                        BandedMatrixQ1<Prec_> A(info.rhs[info.max_level].size() , null.copy(), null.copy() , null.copy(), null.copy() , null.copy(), null.copy(), null.copy(), null.copy(), null.copy());
                        FillMatrix<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(A);
                        DenseVector<Prec_> RHS( info.rhs[info.max_level].size(), Prec_(0.));
                        FillVector<tags::CPU, applications::POISSON, boundary_types::DIRICHLET_NEUMANN>::value(RHS);
                        DenseVector<Prec_> result(n, Prec_(0));

#ifdef HONEI_SSE
                        Multigrid<tags::CPU::SSE, tags::CPU::SSE, NONE, JAC, CYCLE::V, MIXED>::value(A, RHS, result, (unsigned long)11, std::numeric_limits<double>::epsilon(), info);
#else
                        Multigrid<tags::CPU, tags::CPU, NONE, JAC, CYCLE::V, MIXED>::value(A, RHS, result, (unsigned long)11, std::numeric_limits<double>::epsilon(), info);
#endif

                        //write result to scalarfield:
                        unsigned long row_index(0);
                        unsigned long column_index(0);
                        unsigned long index(0);

                        for (; index < info.rhs[info.max_level].size() ; ++index)
                        {
                            (*_u)(row_index, column_index) = result[index] *100;
                            if((index + 1) % _root_n == 0)
                            {
                                ++row_index;
                                column_index = 0;
                            }
                            else
                                ++column_index;
                        }
                    }
                    break;
            }

        }


        void do_timestep(void)
        {
            _update_scenario();
            //do nothing, as in init() the solution is provided
            ++_timestep;
        }

        void render(HONEI_UNUSED bool show_ground, bool use_quads, bool enable_alpha_blending, bool show_water, float alpha)
        {
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
                    glTranslatef(0,0,-200);
                    glBegin(GL_QUADS);
                    for(unsigned int i = 0; i < _root_n-1; ++i)
                    {
                        for(unsigned int j = 0; j <_root_n-1; ++j)
                        {
                            glColor4f(0.0, 0.0, 1.0, alpha);
                            glVertex3d(i,j,((*_u)[j][i]));
                            glColor4f(0.0, 1.0, 1.0, alpha);
                            glVertex3d(i+1,j,((*_u)[j][i+1]));
                            glVertex3d(i+1,j+1,((*_u)[j+1][i+1]));
                            glVertex3d(i,j+1,((*_u)[j+1][i]));
                        }
                    }
                    glEnd();
                }
                else
                {
                    glBegin(GL_TRIANGLE_STRIP);
                    for(unsigned int i = 0; i <  _root_n-1; ++i)
                    {
                        for(unsigned int j = 0; j <  _root_n; j++)
                        {
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i,j, ((*_u)[j][i]));
                            glColor4f(0.0, 0.0, 1.0,  alpha);
                            glVertex3d(i+1,j,((*_u)[j][i+1]));
                        }
                        ++i;
                        if (i >=  _root_n-1)
                            break;
                        for(int j2 =  _root_n-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, ((*_u)[j2][i]));
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i+1,j2, ((*_u)[j2][i+1]));
                        }
                    }
                    glEnd();
                }
            }
        }
};
#endif
