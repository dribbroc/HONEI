/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/mesh.hh>
#include <iostream>
#include <deque>
#include <list>
#include <queue>


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshTest:
    public BaseTest
{
    public:
        MeshTest(const std::string & tag) :
            BaseTest("MeshTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            fem::Mesh<> m(3);

            fem::Mesh<fem::Topology<IndexType_, OT_, IT_> > m2(3);
        }
};
MeshTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > topology_test_cpu_v_v("std::vector, std::vector");
MeshTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > topology_test_cpu_d_v("std::deque, std::vector");
MeshTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > topology_test_cpu_v_d("std::vector, std::deque");
MeshTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > topology_test_cpu_d_d("std::deque, std::deque");
MeshTest<tags::CPU, unsigned long, std::list, std::vector<unsigned long> > topology_test_cpu_l_v("std::list, std::vector");
MeshTest<tags::CPU, unsigned long, std::vector, std::list<unsigned long> > topology_test_cpu_v_l("std::vector, std::list");
MeshTest<tags::CPU, unsigned long, std::list, std::list<unsigned long> > topology_test_cpu_l_l("std::list, std::list");
MeshTest<tags::CPU, unsigned long, std::deque, std::list<unsigned long> > topology_test_cpu_d_l("std::deque, std::list");
MeshTest<tags::CPU, unsigned long, std::list, std::deque<unsigned long> > topology_test_cpu_l_d("std::list, std::deque");
MeshTest<tags::CPU, unsigned long, std::queue, std::queue<unsigned long> > topology_test_cpu_q_q("std::queue, std::queue");
MeshTest<tags::CPU, unsigned long, std::deque, std::queue<unsigned long> > topology_test_cpu_d_q("std::deque, std::queue");
MeshTest<tags::CPU, unsigned long, std::queue, std::deque<unsigned long> > topology_test_cpu_q_d("std::queue, std::deque");
MeshTest<tags::CPU, unsigned long, std::list, std::queue<unsigned long> > topology_test_cpu_l_q("std::list, std::queue");
MeshTest<tags::CPU, unsigned long, std::queue, std::list<unsigned long> > topology_test_cpu_q_l("std::queue, std::list");
MeshTest<tags::CPU, unsigned long, std::queue, std::vector<unsigned long> > topology_test_cpu_q_v("std::queue, std::vector");
MeshTest<tags::CPU, unsigned long, std::vector, std::queue<unsigned long> > topology_test_cpu_v_q("std::vector, std::queue");
