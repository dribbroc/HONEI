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
#include <honei/mpi/matrix_io_mpi.hh>
#include <honei/mpi/vector_io_mpi.hh>
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
        unsigned long _min;

    public:
        MGSolverTest(const std::string & tag, std::string filebase, unsigned long min) :
            BaseTest("MGSolverTest<" + tag + ">")
        {
            register_tag(Tag_::name);
            _file = filebase;
            _min = min;
        }

        virtual void run() const
        {
            Configuration::instance()->set_value("ell::threads", 2);
            Configuration::instance()->set_value("mpi::min_part_size", _min);

            unsigned long max_level(8);
            unsigned long min_level(1);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/";
            file += _file;

            MGData<SparseMatrixELLMPI<double>, DenseVectorMPI<double>, SparseMatrixELLMPI<double>, DenseVectorMPI<double>, double >  data_mpi(MGUtil<Tag_,
                                                                                            SparseMatrixELLMPI<double>,
                                                                                            DenseVectorMPI<double>,
                                                                                            SparseMatrixELLMPI<double>,
                                                                                            DenseVectorMPI<double>,
                                                                                            MatrixIOMPI_ELL<io_formats::ELL>,
                                                                                            VectorIOMPI<io_formats::EXP>,
                                                                                            double>::load_data(file, max_level, double(0.7), "jac"));


            MGUtil<Tag_,
                SparseMatrixELLMPI<double>,
                DenseVectorMPI<double>,
                SparseMatrixELLMPI<double>,
                DenseVectorMPI<double>,
                MatrixIOMPI_ELL<io_formats::ELL>,
                VectorIOMPI<io_formats::EXP>,
                double>::configure(data_mpi, 100, 100, 4, 4, min_level, double(1e-8));

            OperatorList ol_mpi(
                    MGCycleCreation<Tag_,
                    methods::CYCLE::V::STATIC,
                    BiCGStabSolver<Tag_, methods::VAR>,
                    BiCGStabSmoother<Tag_>,
                    //RISmoother<Tag_>,
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
                std::cout<<"min part size: "<<_min<<std::endl;
                std::cout<<"processes: "<<mpi::mpi_comm_size()<<" "<<min_level<<" "<<max_level<<std::endl;
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

            data_mpi.x.at(MGDataIndex::internal_index_A(max_level)).lock(lm_read_only);
            ref.lock(lm_read_only);
            for(unsigned long i(0) ; i < ref.local_size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(data_mpi.x.at(MGDataIndex::internal_index_A(max_level))[i], ref[i], eps);
            data_mpi.x.at(MGDataIndex::internal_index_A(max_level)).unlock(lm_read_only);
            ref.unlock(lm_read_only);

            data_mpi.A.clear();
            data_mpi.P.clear();
            data_mpi.b.clear();
            data_mpi.x.clear();
            data_mpi.c.clear();
            data_mpi.d.clear();
            data_mpi.temp_0.clear();
            data_mpi.temp_1.clear();
            data_mpi.prolmat.clear();
            data_mpi.resmat.clear();
        }
};
MGSolverTest<tags::CPU::Generic> generic_mg_solver_test("double", "poisson_advanced2/q2_sort_0/", 1);
#ifdef HONEI_SSE
MGSolverTest<tags::CPU::SSE> mg0_solver_test_cpu("double", "poisson_advanced2/q2_sort_2/", 1);
#else
MGSolverTest<tags::CPU> mg_solver_test_cpu("double", "poisson_advanced2/sort_2/", 1);
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::GPU::CUDA> cuda_mg_solver_test("double", "poisson_advanced2/q2_sort_0/", 1);
#endif
#endif
