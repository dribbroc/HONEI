/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/math/operator.hh>
#include <honei/math/restriction.hh>
#include <honei/math/prolongation.hh>
#include <honei/math/ir.hh>
#include <honei/math/methods.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename OuterTag_, typename InnerTag_, typename DTOuter_, typename DTInner_>
class IRSolverTest:
    public BaseTest
{
    public:
        IRSolverTest(const std::string & tag) :
            BaseTest("IRSolverTest<" + tag + ">")
        {
            register_tag(OuterTag_::name);
        }

        virtual void run() const
        {
            unsigned long levels(4);
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson_advanced/sort_0/";
            MGData<SparseMatrixELL<DTInner_>, DenseVector<DTInner_>, DenseVector<DTInner_> >  data(MGUtil<InnerTag_,
                                                                                            SparseMatrixELL<DTInner_>,
                                                                                            DenseVector<DTInner_>,
                                                                                            DenseVector<DTInner_>,
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            DTInner_>::load_data(file, levels, DTInner_(0.7), "jac"));
            MGUtil<InnerTag_,
                SparseMatrixELL<DTInner_>,
                DenseVector<DTInner_>,
                DenseVector<DTInner_>,
                io_formats::ELL,
                io_formats::EXP,
                DTInner_>::configure(data, 3, 100, 4, 4, 1, DTInner_(1e-8));

            OperatorList ol(
                    MGCycleCreation<InnerTag_,
                    methods::CYCLE::V::STATIC,
                    CG<InnerTag_, methods::VAR>,
                    RISmoother<InnerTag_>,
                    Restriction<InnerTag_, methods::PROLMAT>,
                    Prolongation<InnerTag_, methods::PROLMAT>,
                    DTInner_>::value(data)
                    );

            std::string a_file(HONEI_SOURCEDIR);
            a_file += "/honei/math/testdata/poisson_advanced/sort_0/A_";
            a_file += stringify(levels);
            a_file += ".ell";
            SparseMatrixELL<DTOuter_> A(MatrixIO<io_formats::ELL>::read_matrix(a_file, DTOuter_(0)));

            std::string b_file(HONEI_SOURCEDIR);
            b_file += "/honei/math/testdata/poisson_advanced/sort_0/rhs_";
            b_file += stringify(levels);
            DenseVector<DTOuter_> b(VectorIO<io_formats::EXP>::read_vector(b_file, DTOuter_(0)));

            std::string x_file(HONEI_SOURCEDIR);
            x_file += "/honei/math/testdata/poisson_advanced/sort_0/init_";
            x_file += stringify(levels);
            DenseVector<DTOuter_> x(VectorIO<io_formats::EXP>::read_vector(x_file, DTOuter_(0)));

            unsigned long used(0);
            IRSolver<OuterTag_, InnerTag_, Norm<vnt_l_two, true, OuterTag_> >::value(A, b, x, data, ol, double(1e-8), 100ul, used);

            std::cout << used << std::endl;
            std::cout << data.used_iters << std::endl;
            std::cout << data.used_iters_coarse << std::endl;

            std::string reffile(HONEI_SOURCEDIR);
            reffile += "/honei/math/testdata/poisson_advanced/sort_0/sol_";
            reffile += stringify(levels);
            DenseVector<double> ref(VectorIO<io_formats::EXP>::read_vector(reffile, double(0)));
            double base_digits(1);
            double additional_digits(2);

            double base_eps(1 / pow(10, base_digits));
            double add_eps(base_eps / pow(10, additional_digits));

            double m((add_eps - base_eps) / double(4));
            double ba(base_eps - (double(4) * m));

            double eps(m * sizeof(double) + ba);
            eps *= double(8);

            x.lock(lm_read_only);
            ref.lock(lm_read_only);
            for(unsigned long i(0) ; i < ref.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(x[i], ref[i], eps);
        }
};
#if defined HONEI_CUDA && defined HONEI_SSE
#ifdef HONEI_CUDA_DOUBLE
IRSolverTest<tags::CPU::MultiCore::SSE, tags::GPU::CUDA, double, float> ir_solver_test_cpu_gpu("ir mcsse/cuda double/float");
#endif
#endif
