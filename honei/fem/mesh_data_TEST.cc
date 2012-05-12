/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/unittest.hh>
#include <honei/fem/mesh_data.hh>
#include <iostream>
#include <deque>
#include <list>


using namespace honei;
using namespace tests;
using namespace std;

template<typename Tag_, typename IndexType_, template<typename, typename> class OT_, typename IT_>
class MeshDataTest:
    public BaseTest
{
    public:
        MeshDataTest(const std::string & tag) :
            BaseTest("MeshDataTest<" + tag + ">")
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

            fem::MeshData<fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > > data(m3);
            DenseVector<double> vertex_attr = data.get_attribute_of_type_1(0);
            TEST_CHECK_EQUAL(vertex_attr.size(), 6u);
        }
};
MeshDataTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > meshdata_test_cpu_v_v("std::vector, std::vector");
MeshDataTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > meshdata_test_cpu_d_v("std::deque, std::vector");
MeshDataTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > meshdata_test_cpu_v_d("std::vector, std::deque");
MeshDataTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > meshdata_test_cpu_d_d("std::deque, std::deque");
