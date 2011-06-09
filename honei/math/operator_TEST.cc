/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/math/operator.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_>
class OperatorTest:
    public BaseTest
{
    public:
        OperatorTest(const std::string & tag) :
            BaseTest("OperatorTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/poisson_advanced/sort_0/A_7.ell";
            SparseMatrixELL<double> system(MatrixIO<io_formats::ELL>::read_matrix(filename, double(0)));

            DenseVector<double> result_1(system.rows());
            DenseVector<double> result_2(system.rows());
            DenseVector<double> rhs(system.rows(), double(1));
            DenseVector<double> x(system.rows(), double(0));

            OperatorList defop;
            defop.push_back(new DefectOperator<Tag_, SparseMatrixELL<double>, DenseVector<double> >(result_1, rhs, system, x));

            defop.value();

            Defect<Tag_>::value(result_2, rhs, system, x);

            for(unsigned long i(0) ; i < x.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(result_2[i], result_1[i], std::numeric_limits<double>::epsilon());

            Operator* sumop;
            sumop = new SumOperator<Tag_, DenseVector<double> >(result_1, rhs);
            Sum<Tag_>::value(result_2, rhs);

            sumop->value();

            for(unsigned long i(0) ; i < x.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(result_2[i], result_1[i], std::numeric_limits<double>::epsilon());

            Operator* axpyop;
            axpyop = new ScaledSumOperator<Tag_, double>(result_1, rhs, double(2));
            ScaledSum<Tag_>::value(result_2, rhs, double(2));

            axpyop->value();

            for(unsigned long i(0) ; i < x.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(result_2[i], result_1[i], std::numeric_limits<double>::epsilon());

            Operator* normop;
            double norm_1;
            normop = new NormOperator<Tag_, vnt_l_two, double, true>(norm_1, rhs);
            double norm_2 = Norm<vnt_l_two, true, Tag_>::value(rhs);

            normop->value();

            TEST_CHECK_EQUAL_WITHIN_EPS(norm_1, norm_2, std::numeric_limits<double>::epsilon());

            Operator* cgop;
            std::string filename2(HONEI_SOURCEDIR);
            filename2 += "/honei/math/testdata/poisson_advanced/sort_0/rhs_7";
            DenseVector<double> b(VectorIO<io_formats::EXP>::read_vector(filename2, double(0)));

            std::string filename3(HONEI_SOURCEDIR);
            filename3 += "/honei/math/testdata/poisson_advanced/sort_0/init_7";
            DenseVector<double> x1(VectorIO<io_formats::EXP>::read_vector(filename3, double(0)));
            DenseVector<double> x2(x1.copy());

            unsigned long used_iters(0);
            cgop = new SolverOperator<CG<Tag_, methods::NONE>, SparseMatrixELL<double>, DenseVector<double> >(system, b, x1, 1000ul, used_iters, double(1e-8));
            cgop->value();

            unsigned long used_iters2(0);
            CG<Tag_, methods::NONE>::value(system, b, x2, 1000ul, used_iters2, double(1e-8));

            TEST_CHECK_EQUAL(used_iters, used_iters2);
            for(unsigned long i(0) ; i < x1.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(x1[i], x2[i], std::numeric_limits<double>::epsilon());
            }

            Operator* riop;

            DenseVector<double> x3(VectorIO<io_formats::EXP>::read_vector(filename3, double(0)));
            DenseVector<double> x4(x3.copy());

            DenseVector<double> t1(x3.size());
            DenseVector<double> t2(x3.size());

            DenseVector<double> diag_inverted(x3.size(), double(0));
            for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
            {
                    diag_inverted[i] = double(0.7)/system(i, i);
            }
            riop = new SmootherOperator<RISmoother<Tag_>, SparseMatrixELL<double>, DenseVector<double>, DenseVector<double> >(system, diag_inverted, b, x3, t1, t2, 1000ul);
            riop->value();

            RISmoother<Tag_>::value(system, diag_inverted, b, x4, t1, t2, 1000ul);

            for(unsigned long i(0) ; i < x4.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(x3[i], x4[i], std::numeric_limits<double>::epsilon());
            }

            Operator* transop;
            string filename_fine(HONEI_SOURCEDIR);
            filename_fine += "/honei/math/testdata/poisson_advanced/sort_0/sol_8";
            string filename_coarse(HONEI_SOURCEDIR);
            filename_coarse += "/honei/math/testdata/poisson_advanced/sort_0/sol_7";

            DenseVector<double> fine_1(VectorIO<io_formats::EXP>::read_vector(filename_fine, double(0)));
            DenseVector<double> fine_2(VectorIO<io_formats::EXP>::read_vector(filename_fine, double(0)));
            DenseVector<double> coarse(VectorIO<io_formats::EXP>::read_vector(filename_coarse, double(0)));

            string filename_prolmat(HONEI_SOURCEDIR);
            filename_prolmat += "/honei/math/testdata/poisson_advanced/sort_0/prol_8.ell";
            SparseMatrixELL<double> prolmat(MatrixIO<io_formats::ELL>::read_matrix(filename_prolmat, double(0)));

            transop = new TransferOperator< Prolongation<Tag_, methods::PROLMAT>, SparseMatrixELL<double>, DenseVector<double> >(fine_1, coarse, prolmat);
            transop->value();

            DenseVector<unsigned long> dummy(1ul);
            Prolongation<Tag_, methods::PROLMAT>::value(fine_2, coarse, dummy, prolmat);
            for(unsigned long i(0) ; i < fine_1.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(fine_1[i], fine_2[i], std::numeric_limits<double>::epsilon());
            }

            delete transop;
            delete riop;
            delete cgop;
            delete axpyop;
            delete sumop;
            delete normop;
        }
};
OperatorTest<tags::CPU> otest_cpu("");
