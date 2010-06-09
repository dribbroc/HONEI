/* vim: set sw=4 sts=4 et foldmethod=syntax : */


#include <honei/math/richardson.hh>
#include <honei/math/preconditioning.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/spai.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class RichardsonTESTArb :
    public BaseTest
{
    public:
        RichardsonTESTArb(const std::string & tag) :
            BaseTest("Richardson solver test (arbitrary -> Jacobi)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            //Read in data:
            std::string dir(HONEI_SOURCEDIR);
            std::string file (dir + "/honei/math/testdata/l2/");
            file += "area51_full_0";
            file += ".ell";
            SparseMatrixELL<DT1_> A(MatrixIO<io_formats::ELL, SparseMatrixELL<double> >::read_matrix(file, DT1_(0)));

            SparseMatrix<DT1_> precon(A);

            std::string rhs_file(dir + "/honei/math/testdata/l2/");
            rhs_file += "area51_rhs_0.dv";
            DenseVector<DT1_> b(VectorIO<io_formats::DV>::read_vector(rhs_file, DT1_(0)));

            //TODO sparsify
            for(typename SparseMatrix<DT1_>::NonZeroElementIterator i(precon.begin_non_zero_elements()) ; i != precon.end_non_zero_elements() ; ++i)
            {
                *i = (i.row() != i.column()) ? DT1_(0) : DT1_(1./ *i);
            }
            SparseMatrixELL<DT1_> C(precon);

            //set up data and info structures:
            DenseVector<DT1_> result(b.size(), DT1_(0));
            DenseVector<DT1_> t0(result.size());
            DenseVector<DT1_> t1(result.size());
            LSData<SparseMatrixELL<DT1_> , SparseMatrixELL<DT1_> , DT1_> data(&A, &C, &b, &result, &t0, &t1);

            unsigned long iters(23);
            LSInfo info(true, true, 1e-08, DT1_(0.7), 1000ul, iters);
            std::cout << "Finished data I/O." << std::endl;
            Richardson<Tag_, Preconditioning<Tag_, methods::NONE> >::value(data, info);

            std::string sol_file(dir + "/honei/math/testdata/l2/");
            sol_file += "area51_sol_0.dv";
            DenseVector<DT1_> sol(VectorIO<io_formats::DV>::read_vector(sol_file, DT1_(0)));

            std::cout << "#Iterations: " << info.iters << std::endl;
            for(unsigned long i(0) ; i < sol.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.result)[i], sol[i], std::numeric_limits<DT1_>::epsilon() * 1e12);
        }
};
RichardsonTESTArb<tags::CPU, double> rt_arb_cpu("double");


template <typename Tag_, typename DT1_>
class RichardsonSPAITEST :
    public BaseTest
{
    public:
        RichardsonSPAITEST(const std::string & tag) :
            BaseTest("Richardson solver test (SPAI)" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            //Read in data:
            std::string dir(HONEI_SOURCEDIR);
            std::string file (dir + "/honei/math/testdata/l2/");
            file += "area51_full_0";
            file += ".ell";
            SparseMatrixELL<DT1_> A(MatrixIO<io_formats::ELL>::read_matrix(file, DT1_(0)));

            SparseMatrix<DT1_> precon(A);
            SparseMatrix<DT1_> spai(SPAI::value(precon, 0.1, 25, 25));

            SparseMatrixELL<DT1_> C(spai);

            std::string rhs_file(dir + "/honei/math/testdata/l2/");
            rhs_file += "area51_rhs_0.dv";
            DenseVector<DT1_> b(VectorIO<io_formats::DV>::read_vector(rhs_file, DT1_(0)));

            //set up data and info structures:
            DenseVector<DT1_> result(b.size(), DT1_(0));
            DenseVector<DT1_> t0(result.size());
            DenseVector<DT1_> t1(result.size());
            LSData<SparseMatrixELL<DT1_> , SparseMatrixELL<DT1_> , DT1_> data(&A, &C, &b, &result, &t0, &t1);

            unsigned long iters(23);
            LSInfo info(true, true, 1e-08, DT1_(0.7), 1000ul, iters);
            std::cout << "Finished data I/O." << std::endl;
            Richardson<Tag_, Preconditioning<Tag_, methods::NONE> >::value(data, info);

            std::string sol_file(dir + "/honei/math/testdata/l2/");
            sol_file += "area51_sol_0.dv";
            DenseVector<DT1_> sol(VectorIO<io_formats::DV>::read_vector(sol_file, DT1_(0)));

            std::cout << "#Iterations: " << info.iters << std::endl;
            for(unsigned long i(0) ; i < sol.size() ; ++i)
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.result)[i], sol[i], std::numeric_limits<DT1_>::epsilon() * 1e12);
        }
};
RichardsonSPAITEST<tags::CPU, double> rt_spai_cpu("double");
