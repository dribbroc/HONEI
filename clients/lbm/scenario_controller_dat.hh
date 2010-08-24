/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef LBM_SCENARIO_CONTROLLER_DAT_HH
#define LBM_SCENARIO_CONTROLLER_DAT_HH

#include <GL/glut.h>
#include <clients/lbm/scenario_controller_base.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

template<typename Tag_, typename Prec_> class ScenarioControllerDat :
    public ScenarioControllerBase
{
    private:
        int scenario_id;

        unsigned long _dwidth, _dheight, _timestep;
        DenseMatrix<Prec_>* _h;

        void _update_scenario()
        {
        }

    public:
        ScenarioControllerDat(int scen_id) :
            scenario_id(scen_id),
            _h(0)
    {
        srand(time(NULL));
    }
        virtual ~ScenarioControllerDat()
        {
            delete _h;
        }

        static int get_precision(int scen_id)
        {
            return 0; // todo return the correct accuracy (0(float) or 1(double))
        }

        void init(void)
        {
            _timestep = 0;
            // read out h and w from out0.dat and alloc matrix
            //
            std::string line;
            std::ifstream file("out0.dat");
            if (!file.is_open())
                throw InternalError("Root file out0.dat not found!");

            getline(file, line);

            std::string buf; // Have a buffer string
            std::stringstream ss(line); // Insert the string into a stream
            std::vector<std::string> tokens; // Create vector to hold our words
            while (ss >> buf)
                tokens.push_back(buf);

            _dheight = atoi(tokens.at(1).c_str());
            _dwidth = atoi(tokens.at(2).c_str());
            std::string title("GnuPlot data " + stringify(_dheight) + " x " + stringify(_dwidth));
            glutSetWindowTitle(title.c_str());
            _h = new DenseMatrix<Prec_>(_dheight, _dwidth, Prec_(0.05));
        }


        void do_timestep(void)
        {
            _update_scenario();
            // read out data from out"timestep".dat
            std::string line;
            std::string filename("out" + stringify(_timestep) + ".dat");
            std::ifstream file(filename.c_str());
            if (!file.is_open())
                _timestep = 0;
            else
            {
                // skip header line
                getline(file, line);
                while (!file.eof())
                {
                    getline(file, line);
                    std::string buf; // Have a buffer string
                    std::stringstream ss(line); // Insert the string into a stream
                    std::vector<std::string> tokens; // Create vector to hold our words
                    while (ss >> buf)
                        tokens.push_back(buf);

                    if (tokens.size() == 3)
                    {
                        unsigned long x(atoi(tokens.at(0).c_str()));
                        unsigned long y(atoi(tokens.at(1).c_str()));
                        Prec_ h(atof(tokens.at(2).c_str()));
                        (*_h)(y, x) = h;
                    }
                }
                ++_timestep;
            }
        }

        void render(bool /*show_ground*/, bool use_quads, bool enable_alpha_blending, bool show_water, float alpha)
        {
            glScalef(1.0f, 1.0f, 100.0f);
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
                            glVertex3d(i,j,((*_h)(j,i) ));
                            glColor4f(0.0, 1.0, 1.0, alpha);
                            glVertex3d(i+1,j,((*_h)[j][i+1] ));
                            glVertex3d(i+1,j+1,((*_h)[j+1][i+1] ));
                            glVertex3d(i,j+1,((*_h)[j+1][i] ));
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
                            glVertex3d(i,j, ((*_h)[j][i]));
                            glColor4f(0.0, 0.0, 1.0,  alpha);
                            glVertex3d(i+1,j,((*_h)[j][i+1] ));
                        }
                        ++i;
                        if (i >=  _dwidth-1)
                            break;
                        for(int j2 =  _dheight-2; j2 >= 0; --j2)
                        {
                            glVertex3d(i,j2, ((*_h)[j2][i]));
                            glColor4f(0.0, 1.0, 1.0,  alpha);
                            glVertex3d(i+1,j2, ((*_h)[j2][i+1] ));
                        }
                    }
                    glEnd();
                }
            }
        }
};
#endif
