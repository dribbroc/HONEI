/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef LBM_GUARD_POST_PROCESSING_HH
#define LBM_GUARD_POST_PROCESSING_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/swe/directions.hh>
#include <iostream>
#include <fstream>

namespace honei
{
    namespace lbm
    {
        struct LBMPostProcessing
        {
            template<typename ResPrec_>
                static inline void value(DenseMatrix<ResPrec_> & height,
                        unsigned long every,
                        unsigned long d_width,
                        unsigned long d_height,
                        ResPrec_ d_x,
                        ResPrec_ d_y,
                        unsigned long solve_time,
                        std::string filename = "")
                {
                    CONTEXT("When processing RelaxSolver post output:");
                    if(solve_time % every == 0)
                    {
                        std::ofstream file;
                        if (filename.length() == 0)
                            filename = "out" + stringify(solve_time) + ".dat";
                        file.open(filename.c_str());

                        ///Write header to file:
                        std::string header = "# " + stringify(d_height) + " " + stringify(d_width) + "\n" + "\n";
                        file << header;

                        ///For every column, traverse elements and write record to file:
                        for(unsigned long x = 0; x < d_width; ++x)
                        {
                            for(unsigned long y = 0; y < d_height; ++y)
                            {
                                std::string record = stringify(x * d_x) + " " + stringify(y * d_y) + " " + stringify((height)(y,x)) + "\n";
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
}
#endif
