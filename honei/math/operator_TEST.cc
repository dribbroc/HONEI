/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <operator.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/matrix_io.hh>

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
            defop.push_back(new DefectOperator<Tag_, SparseMatrixELL<double>, double>(result_1, rhs, system, x));

            defop.value();

            Defect<Tag_>::value(result_2, rhs, system, x);

            for(unsigned long i(0) ; i < x.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS(result_2[i], result_1[i], std::numeric_limits<double>::epsilon());

            Operator* sumop;
            sumop = new SumOperator<Tag_, double>(result_1, rhs);
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

            delete axpyop;
            delete sumop;
            delete normop;
        }
};
OperatorTest<tags::CPU> oest_cpu("cpu double");
