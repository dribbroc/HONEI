/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/math/restriction.hh>
#include <honei/math/prolongation.hh>
#include <honei/math/operator.hh>
#include <honei/math/cg.hh>
#include <honei/math/bicgstab.hh>
#include <honei/math/ri.hh>
#include <honei/math/mg.hh>
#include <honei/math/methods.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>
#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/superlu.hh>
#include <honei/util/time_stamp.hh>

using namespace honei;
using namespace tests;
using namespace std;


template<typename Tag_>
class MGSolverTest:
    public BaseTest
{
    private:
        std::string _file;

    public:
        MGSolverTest(const std::string & tag, std::string filebase) :
            BaseTest("MGSolverTest<" + tag + ">")
        {
            register_tag(Tag_::name);
            _file = filebase;
        }

        virtual void run() const
        {
            unsigned long max_level(6);
            unsigned long min_level(1);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/";
            file += _file;

            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, DenseVector<double> >  data(MGUtil<Tag_,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            double>::load_data(file, max_level, double(0.7), "jac"));

            // create MPI cycle
            std::vector<SparseMatrixELLMPI<double> > A;
            std::vector<SparseMatrixELLMPI<double> > Res;
            std::vector<SparseMatrixELLMPI<double> > Prol;
            std::vector<DenseVectorMPI<double> > P;
            std::vector<DenseVectorMPI<double> > b;
            std::vector<DenseVectorMPI<double> > x;
            std::vector<DenseVectorMPI<double> > c;
            std::vector<DenseVectorMPI<double> > d;
            std::vector<DenseVectorMPI<double> > t0;
            std::vector<DenseVectorMPI<double> > t1;

            for(unsigned long i(0); i < data.A.size(); ++i)
            {
                SparseMatrix<double> as(data.A.at(i));
                SparseMatrixELLMPI<double> am(as);
                A.push_back(am);

                DenseVectorMPI<double> a4m(data.P.at(i));
                P.push_back(a4m);
                DenseVectorMPI<double> a5m(data.b.at(i));
                b.push_back(a5m);
                DenseVectorMPI<double> a6m(data.x.at(i));
                x.push_back(a6m);
                DenseVectorMPI<double> a7m(data.c.at(i));
                c.push_back(a7m);
                DenseVectorMPI<double> a8m(data.d.at(i));
                d.push_back(a8m);
                DenseVectorMPI<double> a9m(data.temp_0.at(i));
                t0.push_back(a9m);
                DenseVectorMPI<double> a10m(data.temp_1.at(i));
                t1.push_back(a10m);
            }

            for(unsigned long i(0); i < data.resmat.size(); ++i)
            {
                SparseMatrix<double> a2s(data.resmat.at(i));
                SparseMatrixELLMPI<double> a2m(a2s);
                Res.push_back(a2m);
                SparseMatrix<double> a3s(data.prolmat.at(i));
                SparseMatrixELLMPI<double> a3m(a3s);
                Prol.push_back(a3m);
            }

            for(unsigned long i(0); i < data.prolmat.size(); ++i)
            {
                SparseMatrix<double> a3s(data.prolmat.at(i));
                SparseMatrixELLMPI<double> a3m(a3s);
                Prol.push_back(a3m);
            }

            MGData<SparseMatrixELLMPI<double>, DenseVectorMPI<double>, SparseMatrixELLMPI<double>, DenseVectorMPI<double> > data_mpi(A, Res, Prol, P, b, x, c, d, t0, t1, 1, 1000, 4, 4, 1, 1e-8);


            MGUtil<Tag_,
                SparseMatrixELLMPI<double>,
                DenseVectorMPI<double>,
                SparseMatrixELLMPI<double>,
                DenseVectorMPI<double>,
                io_formats::ELL,
                io_formats::EXP,
                double>::configure(data_mpi, 100, 100, 4, 4, min_level, double(1e-8));

            OperatorList ol_mpi(
                    MGCycleCreation<Tag_,
                    methods::CYCLE::V::STATIC,
                    BiCGStab<Tag_, methods::VAR>,
                    //SuperLU,
                    RISmoother<Tag_>,
                    Restriction<Tag_, methods::PROLMAT>,
                    Prolongation<Tag_, methods::PROLMAT>,
                    double>::value(data_mpi)
                    );



            TimeStamp at, bt;
            at.take();
            MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data_mpi, ol_mpi);
            bt.take();

            if (mpi::mpi_comm_rank() == 0)
            {
                std::cout<<"Used iters: "<<data_mpi.used_iters<<std::endl;
                std::cout<<"Used coarse iters: "<<data_mpi.used_iters_coarse<<std::endl;
                std::cout<<"TOE: "<<bt.total()-at.total()<<std::endl;
            }


            std::string reffile(HONEI_SOURCEDIR);
            reffile += "/honei/math/testdata/";
            reffile += _file;
            reffile += "sol_";
            reffile += stringify(max_level);
            DenseVector<double> refs(VectorIO<io_formats::EXP>::read_vector(reffile, double(0)));
            DenseVectorMPI<double> ref(refs);
            double base_digits(1);
            double additional_digits(2);

            double base_eps(1 / pow(10, base_digits));
            double add_eps(base_eps / pow(10, additional_digits));

            double m((add_eps - base_eps) / double(4));
            double bd(base_eps - (double(4) * m));

            double eps(m * sizeof(double) + bd);
            eps *= double(8);

            data_mpi.x.at(max_level).lock(lm_read_only);
            ref.lock(lm_read_only);
            for(unsigned long i(0) ; i < ref.local_size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(data_mpi.x.at(max_level)[i], ref[i], eps);
            data_mpi.x.at(max_level).unlock(lm_read_only);
            ref.unlock(lm_read_only);
        }
};
MGSolverTest<tags::CPU> mg_solver_test_cpu("double", "poisson_advanced2/sort_0/");
