/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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


#ifndef LIBSWE_GUARD_POST_PROCESSING_HH
#define LIBSWE_GUARD_POST_PROCESSING_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/swe/directions.hh>
#include <iostream>
#include <fstream>

/**
 * \file
 * Implementation of source processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    namespace output_types
    {
        class GNUPLOT
        {
        };
        class HDF5
        {
        };
    }

    template <typename Type__>
    struct PostProcessing
    {
    };

    /** Implementation of a simple postprocessing method generating a file for GNUPLOT>splot
     *  to display the heightfield in one timestep. The fileformat is "default", which means, that
     *  for each row, a block with tripels is written to the file: <x, y, value>. The row-blocks
     *
     *  are separated by empty records within the file. When reading the file without plotting, be careful to
     *  have in mind, that a step of 1 from one row to the other means a step of one in x -direction so we
     *  have the "scalarfield" -point of view and not the "matrix" -point of view.
     *  To simply display the file type
     *      gnuplot
     *      splot "filename"
     *  For more details see
     *      gnuplot
     *      help splot
     *  The filenames are consisting of the string "rsolverfield" and the number of the related timestep. The filesuffix is "dat".
     *
     * \param every In order to not have an overload of files, we do not need, one can specify the number by which _solve_time has to be able to be divided by restlessly - obviously this will slow down the computation, so the common use case will be to set every to _solve_time, so that we render only the data for the last timestep.
     *
     * \param height The scalarfield to be dumped.
     * \param d_width The width of the grid.
     * \param d_height The height of the grid.
     * \param solve_time The current timestep.
     **/

    template <>
    struct PostProcessing<output_types::GNUPLOT>
    {
        template<typename ResPrec_>
        static inline void value(DenseMatrix<ResPrec_> & height, unsigned long every, unsigned long d_width, unsigned long d_height, unsigned long solve_time)
        {
            CONTEXT("When processing RelaxSolver post output:");
            if(solve_time % every == 0 || solve_time == 0)
            {
                std::string filename;
                std::ofstream file;
                filename = "out" + stringify(solve_time) + ".dat";
                file.open(filename.c_str());

                ///For every column, traverse elements and write record to file:
                for(unsigned long x = 0; x < d_width; ++x)
                {
                    for(unsigned long y = 0; y < d_height; ++y)
                    {
                        std::string record = stringify(x) + " " + stringify(y) + " " + stringify((height)(y,x)) + "\n";
                        file << record;
                    }
                    //Create empty record after each row:
                    file << "\n";
                }
                file.close();
            }
#ifdef SOLVER_VERBOSE
            std::cout <<"Finished postprocessing." << std::endl;
#endif

        }
    };
}
#endif
