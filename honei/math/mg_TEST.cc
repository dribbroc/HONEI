/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifdef HONEI_GMP
#include <gmpxx.h>
#endif

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


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_>
class MGCycleCreationTest:
    public BaseTest
{
    public:
        MGCycleCreationTest(const std::string & tag) :
            BaseTest("MGCycleCreationTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::vector<SparseMatrixELL<double> > A;
            std::vector<SparseMatrixELL<double> > Res;
            std::vector<SparseMatrixELL<double> > Prol;
            std::vector<SparseMatrixELL<double> > P;
            std::vector<DenseVector<double> > b;
            std::vector<DenseVector<double> > x;
            std::vector<DenseVector<double> > c;
            std::vector<DenseVector<double> > d;
            std::vector<DenseVector<double> > t0;
            std::vector<DenseVector<double> > t1;
            for(unsigned long i(0); i < 5; ++i)
            {
                SparseMatrix<double> a_t(1,1);
                a_t[0][0] = 1.;
                SparseMatrixELL<double> a(a_t);
                A.push_back(a);
                Res.push_back(a.copy());
                Prol.push_back(a.copy());
                P.push_back(a.copy());

                DenseVector<double> dummy(1 , 1);
                b.push_back(dummy);
                x.push_back(dummy.copy());
                c.push_back(dummy.copy());
                d.push_back(dummy.copy());
                t0.push_back(dummy.copy());
                t1.push_back(dummy.copy());
            }

            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, SparseMatrixELL<double>, double > data(A, Res, Prol, P, b, x, c, d, t0, t1, 1, 1000, 4, 4, 1, 1e-8);

            OperatorList ol(
            MGCycleCreation<Tag_,
                              methods::CYCLE::V::STATIC,
                              CG<Tag_, methods::NONE>,
                              RISmoother<Tag_>,
                              Restriction<Tag_, methods::PROLMAT>,
                              Prolongation<Tag_, methods::PROLMAT>,
                              double>::value(data)
                    );

            for(unsigned long i(0) ; i < ol.size() ; ++i)
            {
                std::cout << ol[i]->to_string() << std::endl;
            }
        }
};
MGCycleCreationTest<tags::CPU> mgcycproctest_cpu("double");

template<typename Tag_>
class Q1MGCycleCreationTest:
    public BaseTest
{
    public:
        Q1MGCycleCreationTest(const std::string & tag) :
            BaseTest("Q1MGCycleCreationTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::vector<BandedMatrixQx<Q1Type, double> > A;
            std::vector<SparseMatrixELL<double> > Res;
            std::vector<SparseMatrixELL<double> > Prol;
            std::vector<BandedMatrixQx<Q1Type, double> > P;
            std::vector<DenseVector<double> > b;
            std::vector<DenseVector<double> > x;
            std::vector<DenseVector<double> > c;
            std::vector<DenseVector<double> > d;
            std::vector<DenseVector<double> > t0;
            std::vector<DenseVector<double> > t1;
            for(unsigned long i(0); i < 5; ++i)
            {
                BandedMatrix<double> a_t(9);
                BandedMatrixQx<Q1Type, double> a(a_t);
                SparseMatrix<double> trans(a);
                SparseMatrixELL<double> transell(trans);
                A.push_back(a);
                Res.push_back(transell.copy());
                Prol.push_back(transell.copy());
                P.push_back(a.copy());

                DenseVector<double> dummy(1 , 1);
                b.push_back(dummy);
                x.push_back(dummy.copy());
                c.push_back(dummy.copy());
                d.push_back(dummy.copy());
                t0.push_back(dummy.copy());
                t1.push_back(dummy.copy());
            }

            MGData<BandedMatrixQx<Q1Type, double>, DenseVector<double>, SparseMatrixELL<double>, BandedMatrixQx<Q1Type, double>, double > data(A, Res, Prol, P, b, x, c, d, t0, t1, 1, 1000, 4, 4, 1, 1e-8);

            OperatorList ol(
            MGCycleCreation<Tag_,
                              methods::CYCLE::V::STATIC,
                              CG<Tag_, methods::NONE>,
                              RISmoother<Tag_>,
                              Restriction<Tag_, methods::PROLMAT>,
                              Prolongation<Tag_, methods::PROLMAT>,
                              double>::value(data)
                    );

            for(unsigned long i(0) ; i < ol.size() ; ++i)
            {
                std::cout << ol[i]->to_string() << std::endl;
            }
        }
};
Q1MGCycleCreationTest<tags::CPU> q1mgcycproctest_cpu("double");

template<typename Tag_>
class MGUtilLoadTest:
    public BaseTest
{
    public:
        MGUtilLoadTest(const std::string & tag) :
            BaseTest("MGUtilLoadTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long levels(4);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson_advanced/sort_0/";
            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, DenseVector<double>, double > data(MGUtil<Tag_,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            double>::load_data(file, levels, double(0.7), "jac"));

            std::cout << "-----------------------------------" << std::endl;
            for(unsigned long i(0) ; i < data.A.size() ; ++i)
            {
                std::cout << "A_" << i << " is a " << data.A.at(i).rows() << " x " << data.A.at(i).columns() << " matrix." << std::endl;
                std::cout << "P_" << i << " is a container of size" << data.P.at(i).size() << "." << std::endl;
            }
            std::cout << "-----------------------------------" << std::endl;
            for(unsigned long i(0) ; i < data.prolmat.size() ; ++i)
            {
                std::cout << "Prol_" << i << " is a " << data.prolmat.at(i).rows() << " x " << data.prolmat.at(i).columns() << " matrix." << std::endl;
                std::cout << "Res_" << i << " is a " << data.resmat.at(i).rows() << " x " << data.resmat.at(i).columns() << " matrix." << std::endl;
            }
            std::cout << "-----------------------------------" << std::endl;
            for(unsigned long i(0) ; i < data.x.size() ; ++i)
            {
                std::cout << "b_" << i << " is a vector of size " << data.b.at(i).size() << std::endl;
                std::cout << "x_" << i << " is a vector of size " << data.x.at(i).size() << std::endl;
                std::cout << "c_" << i << " is a vector of size " << data.c.at(i).size() << std::endl;
                std::cout << "d_" << i << " is a vector of size " << data.d.at(i).size() << std::endl;
                std::cout << "temp0_" << i << " is a vector of size " << data.temp_0.at(i).size() << std::endl;
                std::cout << "temp1_" << i << " is a vector of size " << data.temp_1.at(i).size() << std::endl;
            }
            std::cout << "-----------------------------------" << std::endl;
        }
};
MGUtilLoadTest<tags::CPU> mgutilloadtest_cpu("double");

template<typename Tag_, typename DT_>
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
#ifdef HONEI_GMP
    mpf_set_default_prec(512);
#endif
            unsigned long max_level(4);
            unsigned long min_level(1);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/";
            file += _file;

            MGData<SparseMatrixELL<DT_>, DenseVector<DT_>, SparseMatrixELL<DT_>, SparseMatrixELL<DT_>, DT_ >  data(MGUtil<Tag_,
                                                                                            SparseMatrixELL<DT_>,
                                                                                            DenseVector<DT_>,
                                                                                            SparseMatrixELL<DT_>,
                                                                                            SparseMatrixELL<DT_>,
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            DT_>::load_data(file, max_level, DT_(1.), "spai_grote"));
            MGUtil<Tag_,
                SparseMatrixELL<DT_>,
                DenseVector<DT_>,
                SparseMatrixELL<DT_>,
                SparseMatrixELL<DT_>,
                io_formats::ELL,
                io_formats::EXP,
                DT_>::configure(data, 100, 100, 8, 8, min_level, DT_(1e-8));

            OperatorList ol(
                    MGCycleCreation<Tag_,
                    methods::CYCLE::V::STATIC,
                    BiCGStab<Tag_, methods::VAR>,
                    //SuperLU,
                    RISmoother<Tag_>,
                    Restriction<Tag_, methods::PROLMAT>,
                    Prolongation<Tag_, methods::PROLMAT>,
                    DT_>::value(data)
                    );

            std::cout << data.A.size() << std::endl;
            std::cout << data.P.size() << std::endl;
            std::cout << data.prolmat.size() << std::endl;
            std::cout << data.resmat.size() << std::endl;
            std::cout << data.b.size() << std::endl;
            std::cout << data.x.size() << std::endl;
            std::cout << data.temp_0.size() << std::endl;
            std::cout << data.temp_1.size() << std::endl;

            MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data, ol);

            std::cout << data.used_iters << std::endl;
            std::cout << data.used_iters_coarse << std::endl;

            std::string reffile(HONEI_SOURCEDIR);
            reffile += "/honei/math/testdata/";
            reffile += _file;
            reffile += "sol_";
            reffile += stringify(max_level);
            DenseVector<DT_> ref(VectorIO<io_formats::EXP>::read_vector(reffile, DT_(0)));
            double base_digits(1);
            double additional_digits(2);

            double base_eps(1 / pow(10, base_digits));
            double add_eps(base_eps / pow(10, additional_digits));

            double m((add_eps - base_eps) / double(4));
            double b(base_eps - (double(4) * m));

            double eps(m * sizeof(double) + b);
            eps *= double(8);

            data.x.at(MGDataIndex::internal_index_A(max_level)).lock(lm_read_only);
            ref.lock(lm_read_only);
            for(unsigned long i(0) ; i < ref.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(data.x.at(MGDataIndex::internal_index_A(max_level))[i], ref[i], eps);
            data.x.at(MGDataIndex::internal_index_A(max_level)).unlock(lm_read_only);
            ref.unlock(lm_read_only);

            //print_cycle(ol, max_level, min_level);
        }
};
MGSolverTest<tags::CPU, double> mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
//MGSolverTest<tags::CPU::MultiCore, double> mc_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
MGSolverTest<tags::CPU::Generic, double> generic_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, double> mc_generic_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
#ifdef HONEI_GMP
MGSolverTest<tags::CPU::Generic, mpf_class> generic_mg_solver_test_cpu_mpf_class("mpf_class", "poisson_advanced/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, mpf_class> mc_generic_mg_solver_test_cpu_mpf_class("mpf_class", "poisson_advanced/sort_0/");
#endif
#ifdef HONEI_SSE
MGSolverTest<tags::CPU::SSE, double> sse_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
MGSolverTest<tags::CPU::MultiCore::SSE, double> mcsse_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::GPU::CUDA, double> mg_solver_test_gpu("double", "poisson_advanced/sort_0/");
#endif
#endif
#ifdef HONEI_OPENCL
MGSolverTest<tags::OpenCL::CPU, double> ocl_cpu_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::OpenCL::GPU, double> ocl_gpu_mg_solver_test_cpu("double", "poisson_advanced/sort_0/");
#endif
#endif

MGSolverTest<tags::CPU, double> mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
MGSolverTest<tags::CPU::Generic, double> generic_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, double> mc_generic_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
#ifdef HONEI_GMP
MGSolverTest<tags::CPU::Generic, mpf_class> generic_mg_solver_test_cpu_2_mpf_class("mpf_class", "poisson_advanced2/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, mpf_class> mc_generic_mg_solver_test_cpu_2_mpf_class("mpf_class", "poisson_advanced2/sort_0/");
#endif
#ifdef HONEI_SSE
MGSolverTest<tags::CPU::SSE, double> sse_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
MGSolverTest<tags::CPU::MultiCore::SSE, double> mcsse_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::GPU::CUDA, double> mg_solver_test_gpu_2("double", "poisson_advanced2/sort_0/");
#endif
#endif
#ifdef HONEI_OPENCL
MGSolverTest<tags::OpenCL::CPU, double> ocl_cpu_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::OpenCL::GPU, double> ocl_gpu_mg_solver_test_cpu_2("double", "poisson_advanced2/sort_0/");
#endif
#endif

MGSolverTest<tags::CPU, double> mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
MGSolverTest<tags::CPU::Generic, double> generic_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, double> mc_generic_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
#ifdef HONEI_GMP
MGSolverTest<tags::CPU::Generic, mpf_class> generic_mg_solver_test_cpu_4_mpf_class("mpf_class", "poisson_advanced4/sort_0/");
MGSolverTest<tags::CPU::MultiCore::Generic, mpf_class> mc_generic_mg_solver_test_cpu_4_mpf_class("mpf_class", "poisson_advanced4/sort_0/");
#endif
#ifdef HONEI_SSE
MGSolverTest<tags::CPU::SSE, double> sse_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
MGSolverTest<tags::CPU::MultiCore::SSE, double> mcsse_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
#endif
#ifdef HONEI_CUDA
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::GPU::CUDA, double> mg_solver_test_gpu_4("double", "poisson_advanced4/sort_0/");
#endif
#endif
#ifdef HONEI_OPENCL
MGSolverTest<tags::OpenCL::CPU, double> ocl_cpu_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
#ifdef HONEI_CUDA_DOUBLE
MGSolverTest<tags::OpenCL::GPU, double> ocl_gpu_mg_solver_test_cpu_4("double", "poisson_advanced4/sort_0/");
#endif
#endif

template<typename Tag_>
class MGSolverTestQuad:
    public BaseTest
{
    public:
        MGSolverTestQuad(const std::string & tag) :
            BaseTest("MGSolverTestQuad<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long max_level(6);
            unsigned long min_level(1);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson/";
            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double>, DenseVector<double>, double >  data(MGUtil<Tag_,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            double>::load_data(file, max_level, double(0.7), "jac"));
            MGUtil<Tag_,
                SparseMatrixELL<double>,
                DenseVector<double>,
                SparseMatrixELL<double>,
                DenseVector<double>,
                io_formats::ELL,
                io_formats::EXP,
                double>::configure(data, 100, 100, 4, 4, min_level, double(1e-8));


            OperatorList ol(
                    MGCycleCreation<Tag_,
                    methods::CYCLE::W::STATIC,
                    CG<Tag_, methods::NONE>,
                    RISmoother<Tag_>,
                    Restriction<Tag_, methods::PROLMAT>,
                    Prolongation<Tag_, methods::PROLMAT>,
                    double>::value(data)
                    );

            MGSolver<Tag_, Norm<vnt_l_two, true, Tag_> >::value(data, ol);

            std::cout << data.used_iters << std::endl;
            std::cout << data.used_iters_coarse << std::endl;

            std::string reffile(HONEI_SOURCEDIR);
            reffile += "/honei/math/testdata/poisson/sol_";
            reffile += stringify(max_level);
            DenseVector<double> ref(VectorIO<io_formats::EXP>::read_vector(reffile, double(0)));
            double base_digits(1);
            double additional_digits(2);

            double base_eps(1 / pow(10, base_digits));
            double add_eps(base_eps / pow(10, additional_digits));

            double m((add_eps - base_eps) / double(4));
            double b(base_eps - (double(4) * m));

            double eps(m * sizeof(double) + b);
            eps *= double(8);


            data.x.at(MGDataIndex::internal_index_A(max_level)).lock(lm_read_only);
            ref.lock(lm_read_only);
            for(unsigned long i(0) ; i < ref.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(data.x.at(MGDataIndex::internal_index_A(max_level))[i], ref[i], eps);
            data.x.at(MGDataIndex::internal_index_A(max_level)).unlock(lm_read_only);
            ref.unlock(lm_read_only);
        }
};
MGSolverTestQuad<tags::CPU> mg_solver_quad_test_cpu("double");
