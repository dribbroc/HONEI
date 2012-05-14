/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/halo.hh>
#include <honei/fem/mesh.hh>
#include <iostream>
#include <deque>
#include <list>


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class HaloTest:
    public BaseTest
{
    public:
        HaloTest(const std::string & tag) :
            BaseTest("HaloTest<" + tag + ">")
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

            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m3;

            //configure attribute
            unsigned my_attribute_index(fem::MeshAttributeRegistration<fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> >, double>::execute(m3, fem::pl_vertex));

            //add vertices
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_polytope(fem::pl_vertex);
            m3.add_attribute_value(my_attribute_index, double(0));
            m3.add_attribute_value(my_attribute_index, double(0.5));
            m3.add_attribute_value(my_attribute_index, double(1));
            m3.add_attribute_value(my_attribute_index, double(0));
            m3.add_attribute_value(my_attribute_index, double(0.5));
            m3.add_attribute_value(my_attribute_index, double(1));

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
            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m4(m3);

            //init simple halo
            fem::Halo<0> h;

            //add connections
            //
            // *--*--*
            // |0 | 1| m3
            // *--*--*
            //  5   6
            //  |   |
            //  0   1
            // *--*--*
            // |0 | 1| m4
            // *--*--*

            h.add_halo_element_pair(5u, 0u);
            h.add_halo_element_pair(6u, 1u);

            TEST_CHECK_EQUAL(h.size(), 2u);
            TEST_CHECK_EQUAL(h.get_element(0u), 5u);
            TEST_CHECK_EQUAL(h.get_element(1u), 6u);
            TEST_CHECK_EQUAL(h.get_element_counterpart(0u), 0u);
            TEST_CHECK_EQUAL(h.get_element_counterpart(1u), 1u);
        }
};
HaloTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > halo_test_cpu_v_v("std::vector, std::vector");
HaloTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > halo_test_cpu_d_v("std::deque, std::vector");
HaloTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > halo_test_cpu_v_d("std::vector, std::deque");
HaloTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > halo_test_cpu_d_d("std::deque, std::deque");
