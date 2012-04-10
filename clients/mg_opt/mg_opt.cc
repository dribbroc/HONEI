/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#include <honei/math/operator.hh>
#include <honei/math/restriction.hh>
#include <honei/math/prolongation.hh>
#include <honei/math/cg.hh>
#include <honei/math/bicgstab.hh>
#include <honei/math/ri.hh>
#include <honei/math/mg.hh>
#include <honei/math/methods.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/superlu.hh>
#include <honei/util/time_stamp.hh>

#ifdef HONEI_MPI
#include <mpi.h>
#endif


using namespace honei;

template <typename Tag_, typename DT_>
class Driver
{
    public:
    Driver(std::string _file)
    {
        if(typeid(Tag_) == typeid(tags::CPU::MultiCore::SSE))
        {
            //Q1
            Configuration::instance()->set_value("ell::threads", 2);
            // Q2 / Q1t
            //Configuration::instance()->set_value("ell::threads", 16);
        }
        else
        {
            //Q1
            Configuration::instance()->set_value("ell::threads", 1);
            Configuration::instance()->set_value("cuda::product_smell_dv_double", 128);
            // Q2 / Q1t
            //Configuration::instance()->set_value("ell::threads", 2);
            //Configuration::instance()->set_value("cuda::product_smell_dv_double", 256);
        }

        unsigned long max_level(8);
        std::string file(HONEI_SOURCEDIR);
        file += "/honei/math/testdata/";
        file += _file;

        MGData<SparseMatrixELL<DT_>, DenseVector<DT_>, SparseMatrixELL<DT_>, DenseVector<DT_>, DT_ >  data_base(MGUtil<Tag_,
                SparseMatrixELL<DT_>,
                DenseVector<DT_>,
                SparseMatrixELL<DT_>,
                DenseVector<DT_>,
                MatrixIO<io_formats::ELL>,
                VectorIO<io_formats::EXP>,
                DT_>::load_data(file, max_level, DT_(1.0), "jac"));

        double best_time(1000);
        double best_iter(1000);
        unsigned long best_min_level(0);
        unsigned long best_smooth_iters(0);
        double best_damping(0);

        for (unsigned long min_level(1) ; min_level <= 3 ; ++min_level)
        {
            std::cout<<"========================= min_level"<<min_level<<"==========================="<<std::endl;
            for (double damping(0.1) ; damping <=1 ; damping += 0.05)
            {
                for (unsigned long smooth_iters(0) ; smooth_iters < 50 ; smooth_iters+=2)
                {
                    MGData<SparseMatrixELL<DT_>, DenseVector<DT_>, SparseMatrixELL<DT_>, DenseVector<DT_>, DT_> data = data_base.copy();
                    for (unsigned long i(0) ; i < data.P.size() ; ++i)
                        Scale<Tag_>::value(data.P.at(i), damping);

                    MGUtil<Tag_,
                        SparseMatrixELL<DT_>,
                        DenseVector<DT_>,
                        SparseMatrixELL<DT_>,
                        DenseVector<DT_>,
                        MatrixIO<io_formats::ELL>,
                        VectorIO<io_formats::EXP>,
                        DT_>::configure(data, 150, 100, smooth_iters, smooth_iters, min_level, DT_(1e-8));

                    OperatorList ol(
                            MGCycleCreation<Tag_,
                            methods::CYCLE::V::STATIC,
                            BiCGStabSolver<Tag_, methods::VAR>,
                            //SuperLU,
                            RISmoother<Tag_>,
                            Restriction<Tag_, methods::PROLMAT>,
                            Prolongation<Tag_, methods::PROLMAT>,
                            DT_>::value(data)
                            );


                    TimeStamp at, bt;
                    at.take();
                    MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data, ol);
                    bt.take();

                    double time(bt.total()-at.total());
                    if (best_time > time && data.used_iters < 140)
                    {
                        best_time = time;
                        best_min_level = min_level;
                        best_smooth_iters = smooth_iters;
                        best_damping = damping;
                        best_iter = data.used_iters;

                        std::cout<<std::endl;
                        std::cout<<"Best time: "<<best_time<<std::endl;
                        std::cout<<"Best iter: "<<best_iter<<std::endl;
                        std::cout<<"Best min level: "<<best_min_level<<std::endl;
                        std::cout<<"Best smooth iters: "<<best_smooth_iters<<std::endl;
                        std::cout<<"Best damping: "<<best_damping<<std::endl;
                    }
                    std::cout<<"."<<std::flush;
                }
            }
        }

        std::cout<<std::endl;
        std::cout<<"Best Time: "<<best_time<<std::endl;
        std::cout<<"Best min level: "<<best_min_level<<std::endl;
        std::cout<<"Best smooth iters: "<<best_smooth_iters<<std::endl;
        std::cout<<"Best damping: "<<best_damping<<std::endl;
    }
};

int main(int argc, char ** argv)
{
#ifdef HONEI_MPI
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
    //Driver<tags::GPU::CUDA, double> a("poisson_advanced2/sort_2/");
#ifdef HONEI_SSE
    Driver<tags::CPU::MultiCore::SSE, double> a("poisson_advanced2/sort_2/");
#endif
#ifdef HONEI_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
