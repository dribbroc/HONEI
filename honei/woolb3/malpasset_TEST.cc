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
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/swe/post_processing.hh>
#include <honei/math/matrix_io.hh>

#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/woolb3/solver_lbm3.hh>

#ifdef DEBUG
#define SOLVER_VERBOSE
#endif

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class MalpassetTest :
    public TaggedTest<Tag_>
{
    public:
        MalpassetTest(const std::string & type) :
            TaggedTest<Tag_>("malpasset_solver_lbm3_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long timesteps(10);
            DataType_ dx(15);
            DataType_ _tau(1.245);
            DataType_ _dt(0.0045);

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

            Grid<D2Q9, DataType_> grid;
            grid.h = new DenseMatrix<DataType_> (h);
            grid.u = new DenseMatrix<DataType_> (u);
            grid.v = new DenseMatrix<DataType_> (v);
            grid.b = new DenseMatrix<DataType_> (b);
            grid.obstacles = new DenseMatrix<bool> (o);

            grid.d_x = dx;
            grid.d_y = dx;
            grid.d_t = _dt;
            grid.tau = _tau;


            Grid3<DataType_, 9> grid3(*grid.obstacles, *grid.h, *grid.b, *grid.u, *grid.v);
            PackedGrid3<DataType_, 9> pgrid3(grid3);
            SolverLBM3<Tag_, DataType_, 9, lbm::lbm_source_schemes::BED_FULL> solver3(grid3, pgrid3, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            solver3.do_preprocessing();

            for (unsigned long i(0) ; i < timesteps ; ++i)
            {
                solver3.solve();
            }


            //solver3.do_postprocessing();

            //DenseMatrix<DataType_> h_3(grid.h->rows(), grid.h->columns(), DataType_(0));
            grid3.fill_h(*grid.h, *pgrid3.h);

            std::string filename(file + "out_");
            filename += stringify(_dt);
            filename += "_";
            filename += stringify(_tau);
            filename +=".dat";
            PostProcessing<GNUPLOT>::value(*grid.h, 1, h.columns(), h.rows(), 0, filename);

            grid.destroy();
        }

};
//MalpassetTest<tags::CPU, float> malpasset_solver_test_float("float");
//MalpassetTest<tags::CPU, double> malpasset_solver_test_double("double");
