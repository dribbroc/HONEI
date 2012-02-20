/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
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

#pragma once
#ifndef MATH_GUARD_LS_CONFIG_HH
#define MATH_GUARD_LS_CONFIG_HH 1

#include <honei/util/exception.hh>
#include <string>
#include <fstream>
#include <iostream>

namespace honei
{
    namespace ls_config
    {
        /// \todo replace double by DT_
        class SolverObject
        {
            public:
            std::string name;
            std::string type;
        };

        class Smoother : public SolverObject
        {
            public:
            unsigned long iters;
            std::string precision;
        };

        class Solver : public SolverObject
        {
            public:
            unsigned long max_iters;
            double convergence;
            std::string precision;
        };

        class MGSmoother : public Smoother
        {
            public:
            SolverObject coarse;
            SolverObject smoother;
        };

        class SGSmoother : public Smoother
        {
            public:
            SolverObject precon;
            double damping;
        };

        class MGSolver : public Solver
        {
            public:
            SolverObject coarse;
            SolverObject smoother;
        };

        class SGSolver : public Solver
        {
            public:
            SolverObject precon;
            double damping;
        };

        class LSConfig
        {
            public:
            static SolverObject parse_file(std::string filename)
            {
                std::fstream file(filename.c_str(), std::ios_base::in);
                SolverObject solver = create_object(file, "solver");
                file.close();

                return solver;
            }

            static SolverObject create_object(std::fstream & file, std::string object_name)
            {
                std::string line;

                // search for the object section position
                file.seekg(0, std::ios::beg);
                bool found = false;
                while (std::getline(file, line))
                {
                    if (line.empty())
                        continue;

                    if ('#' == line.at(0))
                        continue;

                    if (line.find("[" + object_name + "]") != std::string::npos)
                    {
                        found = true;
                        break;
                    }
                }
                if (found == false)
                    throw InternalError("Error: could not find symbol " + object_name + " in solver config file!");

                std::streampos section_pos(file.tellg());

                SolverObject result;
                std::string type(get_value(file, section_pos, "type"));
                std::cout<<"creating object: "<<object_name<<" with type: "<<type<<std::endl;

                if (type == "mg")
                {
                    MGSolver solver;
                    std::string coarse_type(get_value(file, section_pos, "coarse"));
                    solver.coarse = create_object(file, coarse_type);
                    std::string smoother_type(get_value(file, section_pos, "smoother"));
                    solver.smoother = create_object(file, smoother_type);

                    result = solver;
                }
                else if (type == "cg" || type == "ri")
                {
                    SGSolver solver;
                    std::string precon_type(get_value(file, section_pos, "precon"));
                    if (precon_type == "none" || precon_type == "jac" || precon_type == "spai" || precon_type == "sainv")
                    {
                        // TODO load preconditioner
                    }
                    else
                        solver.precon = create_object(file, precon_type);

                    result = solver;
                }
                else if (type == "mg-smoother")
                {
                    MGSmoother smoother;
                    std::string coarse_type(get_value(file, section_pos, "coarse"));
                    smoother.coarse = create_object(file, coarse_type);
                    std::string smoother_type(get_value(file, section_pos, "smoother"));
                    smoother.smoother = create_object(file, smoother_type);

                    result = smoother;
                }
                else if (type == "ri-smoother" || type == "cg-smoother")
                {
                    SGSmoother smoother;
                    std::string precon_type(get_value(file, section_pos, "precon"));
                    if (precon_type == "none" || precon_type == "jac" || precon_type == "spai" || precon_type == "sainv")
                    {
                        // TODO load preconditioner
                    }
                    else
                        smoother.precon = create_object(file, precon_type);

                    result = smoother;
                }
                else
                    throw InternalError("Error: Type " + type + " not known!");


                result.name = object_name;
                result.type = type;

                return result;
            }

            static std::string get_value(std::fstream & file, std::streampos pos, std::string key)
            {
                std::string line;
                file.seekg(pos, std::ios::beg);
                bool found = false;
                while (std::getline(file, line))
                {
                    if (line.empty())
                        continue;

                    if ('#' == line.at(0))
                        continue;

                    if (line.find("[") != std::string::npos)
                        break;

                    if (line.find(key) != std::string::npos)
                    {
                        found = true;
                        break;
                    }
                }
                if (found == false)
                    throw InternalError("Error: could not find symbol " + key + " in corresponding section!");

                size_t value_begin(line.rfind("="));
                if (value_begin == std::string::npos)
                    throw InternalError("Syntax Error: " + line + "!");
                ++value_begin;

                std::string value(line.substr(value_begin, std::string::npos));

                size_t whitespace(value.find(" "));
                while (whitespace != std::string::npos)
                {
                    value.erase(whitespace, 1);
                    whitespace = value.find(" ");
                }

                return value;
            }
        };
    }
}
#endif
