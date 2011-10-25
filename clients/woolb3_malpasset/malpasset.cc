/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/swe/post_processing.hh>
#include <honei/math/matrix_io.hh>

#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/woolb3/solver_lbm3.hh>


using namespace honei;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

typedef double DataType_;


int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'malpasset dt-value dat-filename'"<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string filename(argv[2]);

    unsigned long timesteps(5000);
    DataType_ dx(15);
    DataType_ _tau(0.51);
    //DataType_ _dt(0.0015);
    DataType_ _dt(atof(argv[1]));

    std::string file(HONEI_SOURCEDIR);
    file += "/honei/woolb3/malpasset/";
    SparseMatrixELL<DataType_> pre_h_b(MatrixIO<io_formats::ELL>::read_matrix(file + "malp_initial.ell", DataType_(0)));
    DenseMatrix<DataType_> h_b(pre_h_b.rows(), pre_h_b.columns());

    SparseMatrixELL<DataType_> pre_pre_s(MatrixIO<io_formats::ELL>::read_matrix(file + "malp_domain.ell", DataType_(0)));
    DenseMatrix<DataType_> pre_s(pre_h_b.rows(), pre_h_b.columns());

    SparseMatrixELL<DataType_> pre_b(MatrixIO<io_formats::ELL>::read_matrix(file + "malp_bottom.ell", DataType_(0)));
    DenseMatrix<DataType_> b(pre_h_b.rows(), pre_h_b.columns());

    unsigned long height(pre_h_b.rows());
    unsigned long width(pre_h_b.columns());
    DenseMatrix<DataType_> h((unsigned long)height, (unsigned long)width, DataType_(0));
    DenseMatrix<DataType_> u((unsigned long)height, (unsigned long)width, DataType_(0));
    DenseMatrix<DataType_> v((unsigned long)height, (unsigned long)width, DataType_(0));
    DenseMatrix<bool> o((unsigned long)height, (unsigned long)width, bool(0));
    DenseMatrix<bool> s((unsigned long)height, (unsigned long)width, bool(0));

    for (unsigned long i = 0; i < height; ++i)
    {
        for (unsigned long j = 0; j < width; ++j)
        {
            h_b[i][j] = pre_h_b(i , j);
            pre_s[i][j] = pre_pre_s(i , j);
            s[i][j] = pre_s[i][j] == 110 ? true : false;
            o[i][j] = pre_s[i][j] == 110 ? true : false;
            b[i][j] = s(i , j) ? float(0) : pre_s(i , j);
            h[i][j] = h_b[i][j] - b[i][j];
        }
    }

    DenseMatrix<DataType_> twenty((unsigned long)height, (unsigned long)width, DataType_(20));
    Sum<tags::CPU>::value(b ,twenty);
    Sum<tags::CPU>::value(h ,twenty);

    Scale<tags::CPU>::value(b, DataType_(1./1000));
    Scale<tags::CPU>::value(h, DataType_(1./1000));
    dx /= 1000;


    Grid3<DataType_, 9> grid3(o, h, b, u, v);




    for ( ; _dt <= 0.15 ; _dt+=0.01 )
    {
        _tau = DataType_(0.51);
        for( ; _tau <= 2 ; _tau+=0.01 )
        {
            PackedGrid3<DataType_, 9> pgrid3(grid3);
            SolverLBM3<tags::CPU, DataType_, 9, lbm::lbm_source_schemes::BED_FULL> solver3(grid3, pgrid3, dx, dx, _dt, _tau);
            solver3.do_preprocessing();

            unsigned long i(0);
            //DenseVector<DataType_> old(pgrid3.h->copy());
            for (bool broken(false) ; (i < timesteps) && (!broken) ; ++i)
            {
                solver3.solve();
                for (unsigned long j(0) ; j < pgrid3.h->size() ; ++j)
                {
                    if ( (*pgrid3.h)[j] > DataType_(0.15))
                    {
                        broken = true;
                        break;
                    }
                    else
                    {
                        //old = pgrid3.h->copy();
                    }
                }
            }

            //if (i > best_iters)
            {
                //std::cout<<std::endl<<"iters: "<<i<<" dt: "<<_dt<<" tau: "<<_tau<<std::endl;
                //best_iters = i;
                /*DenseMatrix<DataType_> output(h.copy());
                  grid3.fill_h(output, old);

                  std::string filename(file + "out_");
                  filename += stringify(i);
                  filename += "-";
                  filename += stringify(_dt);
                  filename += "_";
                  filename += stringify(_tau);
                  filename +=".dat";
                  PostProcessing<GNUPLOT>::value(output, 1, output.columns(), output.rows(), 0, filename);*/


                std::ofstream file;
                //std::string filename(HONEI_SOURCEDIR);
                //filename += "/honei/woolb3/malpasset/";
                //filename += "heatmap.dat";
                file.open(filename.c_str(), std::ios::ate | std::ios::out | std::ios::app);
                //std::string header = "# " + stringify(d_height) + " " + stringify(d_width) + "\n" + "\n";
                //file << header;
                std::string record = stringify(_dt) + " " + stringify(_tau) + " " + stringify(i) + "\n";
                file << record;
                file.close();
            }
            /*else
              {
              std::cout<<i<<" ";
              std::cout.flush();
              }*/
        }
    }
}

