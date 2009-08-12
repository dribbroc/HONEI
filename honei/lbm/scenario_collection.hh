/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LBM_GUARD_SCENARIO_COLLECTION_HH
#define LBM_GUARD_SCENARIO_COLLECTION_HH 1

/**
 * \file
 * Implementation of a central scenario database for use in
 * test/benchmarks/clients within the HONEI LBM framework.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/grid.hh>
#include <honei/swe/volume.hh>
#include <cmath>

using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

class ScenarioCollection
{
    private:
        static const unsigned long _scenario_count = 9;
        static const unsigned long _stable_scenario_count = 6;

    public:
        static std::string get_scenario_descr(unsigned long id)
        {
            switch(id)
            {
                case 0:
                    return "Circular dam break, uncritical";
                    break;

                case 1:
                    return "Circular dam break, uncritical, with cuboidal obstacles";
                    break;

                case 2:
                    return "Partial cuboidal dam break, uncritical";
                    break;

                case 3:
                    return "Circular dam break over uneven bed topography, uncritical";
                    break;

                case 4:
                    return "Circular dam break over uneven bed topography and cylindrical obstacle, uncritical";
                    break;

                case 5:
                    return "Circular dam break over uneven bed topography, CRITICAL(dry states)";
                    break;

                case 6:
                    return "STEADY STATE";
                    break;

                case 7:
                    return "Partial cuboidal dam break, CRITICAL(dry states)";
                    break;

                case 8:
                    return "DRIVEN CAVITY";
                    break;

                default:
                    return "Unspecified or unstable  scenario";
                    break;

            }
        }
        template <typename GridType_, typename DataType_>
        static void get_scenario(unsigned long scen_id, unsigned long grid_height, unsigned long grid_width, Grid<GridType_, DataType_> & target_grid)
        {
            switch (scen_id)
            {
                case 0:
                    {
                        target_grid.description = "Circular dam break, uncritical";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n");
                        target_grid.long_description.append("recommended minimum size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.03), 25, 25);
                        c1.value();
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                    }
                    break;
                case 1:
                    {
                        target_grid.description = "Circular dam break, uncritical, with cuboidal obstacles";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");
                        target_grid.long_description.append("recommended minimum size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.06), 25, 25);
                        c1.value();
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);
                        Cuboid<bool> q2(*target_grid.obstacles, 15, 5, 1, 10, 0);
                        q2.value();
                        Cuboid<bool> q3(*target_grid.obstacles, 40, 5, 1, 10, 30);
                        q3.value();

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                    }
                    break;
                case 2:
                    {
                        target_grid.description = "Partial cuboidal dam break, uncritical";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Scenario constraints:\n");
                        target_grid.long_description.append("size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h_upper = h_upper + b = 0.05\n");
                        target_grid.long_description.append("h_lower = h_lower + b = 0.085\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        unsigned long constraint_g_h = 50;
                        unsigned long constraint_g_w = 50;

                        target_grid.h = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.05));
                        Cuboid<DataType_> reservoir(*target_grid.h, 50, 18, DataType_(0.035), 32, 0);
                        reservoir.value();
                        target_grid.u = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(constraint_g_h, constraint_g_w, false);
                        Cuboid<bool> q2(*target_grid.obstacles, 15, 2, 1, 30, 0);
                        q2.value();
                        Cuboid<bool> q3(*target_grid.obstacles, 40, 2, 1, 30, 20);
                        q3.value();

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;
                    }
                    break;
                case 3:
                    {
                        target_grid.description = "Circular dam break over uneven bed topography, uncritical";
                        target_grid.long_description = target_grid.description;

                        target_grid.long_description.append("\nScenario constraints:\n");
                        target_grid.long_description.append("size: 100x200\n");
                        target_grid.long_description.append("forces: BED_SLOPE, BED_FRICTION\n\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = predefined 2D function f with f_max = 0.04\n");
                        target_grid.long_description.append("h = h - f\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        unsigned long constraint_g_h = 100;
                        unsigned long constraint_g_w = 200;

                        target_grid.h = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.03), 35 ,16);
                        c1.value();

                        target_grid.u = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(constraint_g_h, constraint_g_w, false);

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                        //build up the hill:
                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                double x(j * target_grid.d_x);
                                double y(i * target_grid.d_y);
                                //if(sqrt(y * y + x * x) >= 0.4)
                                    (*target_grid.b)(i , j) = 0.04 * exp((-5.) * (x - 1.) * (x - 1.) - 50. * (y - 0.5) * (y - 0.5));
                            }
                        }

                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                    (*target_grid.h)(i , j) -= (*target_grid.b)(i , j);
                            }
                        }


                    }
                    break;
                case 4:
                    {
                        target_grid.description = "Circular dam break over uneven bed topography and cylindrical obstacle, uncritical";
                        target_grid.long_description = target_grid.description;

                        target_grid.long_description.append("\nScenario constraints:\n");
                        target_grid.long_description.append("size: 100x200\n");
                        target_grid.long_description.append("forces: BED_SLOPE, BED_FRICTION\n\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = predefined 2D function f with f_max = 0.04\n");
                        target_grid.long_description.append("h = h - f\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        unsigned long constraint_g_h = 100;
                        unsigned long constraint_g_w = 100;

                        target_grid.h = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.03), 35 ,16);
                        c1.value();

                        target_grid.u = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(constraint_g_h, constraint_g_w, false);
                        Cylinder<bool> c2(*target_grid.obstacles, true, 10 ,10);
                        c2.value();

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                        //build up the hill:
                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                double x(j * target_grid.d_x);
                                double y(i * target_grid.d_y);
                                //if(sqrt(y * y + x * x) >= 0.4)
                                    (*target_grid.b)(i , j) = 0.04 * exp((-5.) * (x - 1.) * (x - 1.) - 50. * (y - 0.5) * (y - 0.5));
                            }
                        }

                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                    (*target_grid.h)(i , j) -= (*target_grid.b)(i , j);
                            }
                        }


                    }
                    break;
                case 5:
                    {
                        target_grid.description = "Circular dam break over uneven bed topography, CRITICAL(dry states)";
                        target_grid.long_description = target_grid.description;

                        target_grid.long_description.append("\nScenario constraints:\n");
                        target_grid.long_description.append("size: 100x200\n");
                        target_grid.long_description.append("forces: BED_SLOPE, BED_FRICTION\n\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = predefined 2D function f with f_max = 0.08\n");
                        target_grid.long_description.append("h = h - f\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        unsigned long constraint_g_h = 100;
                        unsigned long constraint_g_w = 200;

                        target_grid.h = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.05));
                        Cylinder<DataType_> c1(*target_grid.h, DataType_(0.03), 35 ,16);
                        c1.value();

                        target_grid.u = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(constraint_g_h, constraint_g_w, false);
                        Cylinder<bool> c2(*target_grid.obstacles, true, 10 ,10);
                        c2.value();

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                        //build up the hill:
                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                double x(j * target_grid.d_x);
                                double y(i * target_grid.d_y);
                                //if(sqrt(y * y + x * x) >= 0.4)
                                    (*target_grid.b)(i , j) = 0.08 * exp((-5.) * (x - 1.) * (x - 1.) - 50. * (y - 0.5) * (y - 0.5));
                            }
                        }

                        for(unsigned long i(0) ; i < constraint_g_h ; ++i)
                        {
                            for(unsigned long j(0) ; j < constraint_g_w ; ++j)
                            {
                                    (*target_grid.h)(i , j) -= (*target_grid.b)(i , j);
                            }
                        }


                    }
                    break;
                case 6:
                    {
                        target_grid.description = "STEADY STATE";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n");
                        target_grid.long_description.append("recommended minimum size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);

                        target_grid.d_x = 0.1;
                        target_grid.d_y = 0.1;
                        target_grid.d_t = 0.1;
                        target_grid.tau = 1.;

                    }
                    break;
                case 7:
                    {
                        target_grid.description = "Partial cuboidal dam break, CRITICAL(dry states)";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n\n");

                        target_grid.long_description.append("Scenario constraints:\n");
                        target_grid.long_description.append("size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h_upper = h_upper + b = 0.0\n");
                        target_grid.long_description.append("h_lower = h_lower + b = 0.03\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        unsigned long constraint_g_h = 50;
                        unsigned long constraint_g_w = 50;

                        target_grid.h = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.0));
                        Cuboid<DataType_> reservoir(*target_grid.h, 50, 18, DataType_(0.03), 32, 0);
                        reservoir.value();
                        target_grid.u = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(constraint_g_h, constraint_g_w, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(constraint_g_h, constraint_g_w, false);
                        Cuboid<bool> q2(*target_grid.obstacles, 15, 2, 1, 30, 0);
                        q2.value();
                        Cuboid<bool> q3(*target_grid.obstacles, 40, 2, 1, 30, 20);
                        q3.value();

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;
                    }
                    break;
                case 8:
                    {
                        target_grid.description = "DRIVEN CAVITY";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 0.01\n");
                        target_grid.long_description.append("d_t = 0.01\n");
                        target_grid.long_description.append("recommended minimum size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);

                        target_grid.d_x = 0.01;
                        target_grid.d_y = 0.01;
                        target_grid.d_t = 0.01;
                        target_grid.tau = 1.1;

                    }
                    break;
                case 9:
                    {
                        target_grid.description = "STEADY STATE 2";
                        target_grid.long_description = target_grid.description;
                        target_grid.long_description.append("\nDiscretization (proposal):\n");
                        target_grid.long_description.append("d_x = d_y = 1\n");
                        target_grid.long_description.append("d_t = 1\n");
                        target_grid.long_description.append("recommended minimum size: 50x50\n\n");

                        target_grid.long_description.append("Initial conditions:\n");
                        target_grid.long_description.append("u = 0.\n");
                        target_grid.long_description.append("v = 0.\n");
                        target_grid.long_description.append("b = 0.\n");
                        target_grid.long_description.append("h = h + b = 0.05\n\n");

                        target_grid.long_description.append("Lattice Boltzmann model:\n");
                        target_grid.long_description.append("tau = 1.1\n");
                        target_grid.long_description.append("flow = laminar\n");
                        target_grid.long_description.append("lattice_type = D2Q9 square\n");

                        target_grid.h = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.05));
                        target_grid.u = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.v = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.b = new DenseMatrix<DataType_>(grid_height, grid_width, DataType_(0.));
                        target_grid.obstacles = new DenseMatrix<bool>(grid_height, grid_width, false);

                        target_grid.d_x = 1.;
                        target_grid.d_y = 1.;
                        target_grid.d_t = 1.;
                        target_grid.tau = 1.;

                    }
                    break;
            }
            return;
        }

        static unsigned long get_scenario_count()
        {
            return _scenario_count;
        }

        static unsigned long get_stable_scenario_count()
        {
            return _stable_scenario_count;
        }
};


#endif
