/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/halo.hh>
#include <honei/fem/halo_data.hh>
#include <honei/fem/mesh.hh>
#include <honei/fem/communication.hh>
#include <iostream>
#include <deque>
#include <list>


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class CommunicationTest:
    public BaseTest
{
    public:
        CommunicationTest(const std::string & tag) :
            BaseTest("CommunicationTest<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {

            //##################################################################
            //     0  1
            //   0--1--2     *--*--*
            // 2 | 3|  |4    | 0| 1|
            //   3--4--5     *--*--*
            //    5  6

            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m3(0);

            //configure attribute
            unsigned my_attribute_index(fem::MeshAttributeRegistration<fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> >, double>::execute(m3, fem::pl_face));

            //add vertices
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);

            //add edges
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);
            m3.add_polytope(fem::pl_edge);

            //add faces
            m3.add_polytope(fem::pl_face);
            m3.add_polytope(fem::pl_face);
            m3.add_attribute_value(my_attribute_index, double(42.));
            m3.add_attribute_value(my_attribute_index, double(47.));

            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 0, 0); //v->e is set automagically
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 0, 1);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 1, 1);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 1, 2);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 2, 0);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 2, 3);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 3, 1);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 3, 4);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 4, 2);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 4, 5);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 5, 3);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 5, 4);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 6, 4);
            m3.add_adjacency(fem::pl_edge, fem::pl_vertex, 6, 5);

            m3.add_adjacency(fem::pl_face, fem::pl_edge, 0, 0);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 0, 2);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 0, 3);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 0, 5);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 1, 1);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 1, 3);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 1, 4);
            m3.add_adjacency(fem::pl_face, fem::pl_edge, 1, 6);

            //clone mesh
            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m4(1, m3);
            //alter m4's attribute values
            m4.set_attribute_value(my_attribute_index, 0u, 3333.);
            m4.set_attribute_value(my_attribute_index, 1u, 4444.);

            //init simple halo
            fem::Halo<0, fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > > h(m3, 1);

            //add connections
            //
            // *--*--*
            // |0 | 1| m3
            // *--*--*
            //  |   |
            // *--*--*
            // |0 | 1| m4
            // *--*--*

            h.add_halo_element_pair(0u, 0u);
            h.add_halo_element_pair(1u, 1u);

            TEST_CHECK_EQUAL(h.size(), 2u);
            TEST_CHECK_EQUAL(h.get_element(0u), 0u);
            TEST_CHECK_EQUAL(h.get_element(1u), 1u);
            TEST_CHECK_EQUAL(h.get_element_counterpart(0u), 0u);
            TEST_CHECK_EQUAL(h.get_element_counterpart(1u), 1u);

            fem::HaloData<
                fem::Halo<0, fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > > > hd(h);

            //reference to m4 would have been resolved locally
            fem::Communication<0, fem::com_send_receive, double, tags::CPU>::execute(h, 0u, m4, 0u);
            TEST_CHECK_EQUAL(m3.get_attributes_of_type_1().at(my_attribute_index).at(0), 3333.);
            TEST_CHECK_EQUAL(m3.get_attributes_of_type_1().at(my_attribute_index).at(1), 4444.);
            TEST_CHECK_EQUAL(m4.get_attributes_of_type_1().at(my_attribute_index).at(0), 42.);
            TEST_CHECK_EQUAL(m4.get_attributes_of_type_1().at(my_attribute_index).at(1), 47.);

            fem::Communication<0, fem::com_send_receive, double, tags::CPU>::execute(hd, 0u, m4, 0u);
            TEST_CHECK_EQUAL(m3.get_attributes_of_type_1().at(my_attribute_index).at(0), 42.);
            TEST_CHECK_EQUAL(m3.get_attributes_of_type_1().at(my_attribute_index).at(1), 47.);
            TEST_CHECK_EQUAL(m4.get_attributes_of_type_1().at(my_attribute_index).at(0), 3333.);
            TEST_CHECK_EQUAL(m4.get_attributes_of_type_1().at(my_attribute_index).at(1), 4444.);

        }
};
CommunicationTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > halo_test_cpu_v_v("std::vector, std::vector");
CommunicationTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > halo_test_cpu_d_v("std::deque, std::vector");
CommunicationTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > halo_test_cpu_v_d("std::vector, std::deque");
CommunicationTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > halo_test_cpu_d_d("std::deque, std::deque");
