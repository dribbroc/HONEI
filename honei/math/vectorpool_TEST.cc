/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifdef HONEI_GMP
#include <gmpxx.h>
#endif

#include <honei/math/vectorpool.hh>
#include <honei/util/unittest.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/ri.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_>
class VectorpoolTest:
    public BaseTest
{
    public:
        VectorpoolTest(const std::string & tag) :
            BaseTest("VectorpoolTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::vector<DenseVector<double> > stv = create_vectorpool<DenseVector<double> >(10, 100);

            for(unsigned long i(0) ; i < 10 ; ++i)
            {
                TEST_CHECK_EQUAL(stv.at(i).size(), 100ul);
            }

            std::vector<DenseVector<double> > stv2;
            RISmoother<Tag_>::vectorpool(100, stv2);
            for(unsigned long i(0) ; i < 2 ; ++i)
            {
                TEST_CHECK_EQUAL(stv.at(i).size(), 100ul);
            }

            std::vector<std::vector<DenseVector<double> > > data;
            for(unsigned long i(0) ; i < 10 ; ++i)
            {
                std::vector<DenseVector<double> > stv3;
                RISmoother<Tag_>::vectorpool(i*100, stv3);
                data.push_back(stv3);
            }
            for(unsigned long i(0) ; i < 10 ; ++i)
            {
                for(unsigned long j(0) ; j < 2 ; ++j)
                    TEST_CHECK_EQUAL(data.at(i).at(j).size(), i*100ul);
            }
        }
};

VectorpoolTest<tags::CPU> vptest("dv double");
