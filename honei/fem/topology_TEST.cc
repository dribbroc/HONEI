/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/topology.hh>
#include <iostream>
#include <deque>
#include <list>
#include <queue>
#include <dense_data_wrapper.hh>


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class TopologyTest:
    public BaseTest
{
    public:
        TopologyTest(const std::string & tag) :
            BaseTest("TopologyTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            fem::Topology<> t;

            t.push_back();

            TEST_CHECK_EQUAL(t.size(), 1ul);

            t.at(0).push_back(1123);
            t[0].push_back(878);
            TEST_CHECK_EQUAL(t.at(0).at(0), 1123ul);
            TEST_CHECK_EQUAL(t[0][1], 878ul);

            TEST_CHECK_EQUAL(t.size(), 1ul);

            fem::Topology<unsigned long, OT_, IT_> t2;
        }
};
TopologyTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > topology_test_cpu_v_v("std::vector, std::vector");
TopologyTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > topology_test_cpu_d_v("std::deque, std::vector");
TopologyTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > topology_test_cpu_v_d("std::vector, std::deque");
TopologyTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > topology_test_cpu_d_d("std::deque, std::deque");
TopologyTest<tags::CPU, unsigned long, std::list, std::vector<unsigned long> > topology_test_cpu_l_v("std::list, std::vector");
TopologyTest<tags::CPU, unsigned long, std::vector, std::list<unsigned long> > topology_test_cpu_v_l("std::vector, std::list");
TopologyTest<tags::CPU, unsigned long, std::list, std::list<unsigned long> > topology_test_cpu_l_l("std::list, std::list");
TopologyTest<tags::CPU, unsigned long, std::deque, std::list<unsigned long> > topology_test_cpu_d_l("std::deque, std::list");
TopologyTest<tags::CPU, unsigned long, std::list, std::deque<unsigned long> > topology_test_cpu_l_d("std::list, std::deque");

TopologyTest<tags::CPU, unsigned long, std::vector, fem::DenseDataWrapper<15, unsigned long> > topology_test_cpu_v_ddw("std::vector, DV");
TopologyTest<tags::CPU, unsigned long, std::deque, fem::DenseDataWrapper<15, unsigned long> > topology_test_cpu_d_ddw("std::deque, DV");
TopologyTest<tags::CPU, unsigned long, std::list, fem::DenseDataWrapper<15, unsigned long> > topology_test_cpu_l_ddw("std::list, DV");
