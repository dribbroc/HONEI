/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/math/operator.hh>
#include <honei/math/restriction.hh>
#include <honei/math/prolongation.hh>
#include <honei/math/cg.hh>
#include <honei/math/ri.hh>
#include <honei/math/mg.hh>
#include <honei/math/methods.hh>
#include <honei/util/unittest.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_>
class MGCycleProcessingTest:
    public BaseTest
{
    public:
        MGCycleProcessingTest(const std::string & tag) :
            BaseTest("MGCycleProcessingTest<" + tag + ">")
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

            unsigned long used_iters(0);
            MGData<SparseMatrixELL<double>, DenseVector<double>, SparseMatrixELL<double> > data(A, Res, Prol, P, b, x, c, d, t0, t1, 1000, used_iters, 4, 4, 1, 1e-8);

            OperatorList ol(
            MGCycleProcessing<Tag_,
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
MGCycleProcessingTest<tags::CPU> mgcycproctest_cpu("double");

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
            std::string file(HONEI_SOURCEDIR);
            file += "/honei/math/testdata/poisson_advanced/sort_0/";
            MGData<SparseMatrixELL<double>, DenseVector<double>, int> data(MGUtil<Tag_, //TODO: no int
                                                                                            SparseMatrixELL<double>,
                                                                                            DenseVector<double>,
                                                                                            int, //TODO no int
                                                                                            io_formats::ELL,
                                                                                            io_formats::EXP,
                                                                                            double>::load_data(file, 10));
        }
};
MGUtilLoadTest<tags::CPU> mgutilloadtest_cpu("double");
