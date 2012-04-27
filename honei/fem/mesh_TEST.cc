/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/mesh.hh>
#include <iostream>
#include <deque>
#include <list>
#include <queue>
#include <honei/fem/dense_data_wrapper.hh>


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
            fem::Mesh<> m;

            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m2;

            TEST_CHECK_EQUAL(m.get_inter_topology_relation(0), 2u);
            TEST_CHECK_EQUAL(m.get_inter_topology_relation(1), 0u);
            TEST_CHECK_EQUAL(m.get_inter_topology_relation(2), 1u);
            TEST_CHECK_EQUAL(m2.get_inter_topology_relation(0), 2u);
            TEST_CHECK_EQUAL(m2.get_inter_topology_relation(1), 0u);
            TEST_CHECK_EQUAL(m2.get_inter_topology_relation(2), 1u);

            IT_ neighbours;
            neighbours.push_back(234);
            neighbours.push_back(546734);
            neighbours.push_back(2);

            m2.add_element(neighbours);

            IT_ adjacencies;
            adjacencies.push_back(1);

            m2.add_polytope(fem::pl_vertex, adjacencies);

            m2.add_element();
            m2.add_adjacency(fem::pl_face, 0, 1);
        }
};
MeshTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > topology_test_cpu_v_v("std::vector, std::vector");
MeshTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > topology_test_cpu_d_v("std::deque, std::vector");
MeshTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > topology_test_cpu_v_d("std::vector, std::deque");
MeshTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > topology_test_cpu_d_d("std::deque, std::deque");

//MeshTest<tags::CPU, unsigned long, std::vector, fem::DenseDataWrapper<15, unsigned long> > topology_test_cpu_v_ddw("std::vector, DV");
