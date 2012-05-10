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
            //basic tests
            fem::Mesh<> m;

            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m2;

            TEST_CHECK_EQUAL(m.get_num_levels(), 3ul);
            TEST_CHECK_EQUAL(m2.get_num_levels(), 3ul);

            TEST_CHECK_EQUAL(m2.get_downward_index(fem::pl_vertex), -1);
            TEST_CHECK_EQUAL(m2.get_downward_index(fem::pl_edge), fem::ipi_edge_vertex);
            TEST_CHECK_EQUAL(m2.get_downward_index(fem::pl_face), fem::ipi_face_edge);
            TEST_CHECK_EQUAL(m2.get_downward_index(fem::pl_polyhedron), -1);

            TEST_CHECK_EQUAL(m2.get_upward_index(fem::pl_vertex), fem::ipi_vertex_edge);
            TEST_CHECK_EQUAL(m2.get_upward_index(fem::pl_edge), fem::ipi_edge_face);
            TEST_CHECK_EQUAL(m2.get_upward_index(fem::pl_face), -1);
            TEST_CHECK_EQUAL(m2.get_upward_index(fem::pl_polyhedron), -1);

            //##################################################################
            //     0  1
            //   0--1--2     *--*--*
            // 2 | 3|  |4    | 0| 1|
            //   3--4--5     *--*--*
            //    5  6

            fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> > m3;

            //configure attribute
            unsigned my_attribute(fem::MeshAttributeRegistration<fem::Mesh<fem::rnt_2D, fem::Topology<IndexType_, OT_, IT_> >, double>::execute(m3, fem::pl_vertex));
            std::cout << my_attribute << std::endl;

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

            //testing face-edge access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_1(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_edge, 0));
            TEST_CHECK_EQUAL(test_1.size(), 4ul);
            TEST_CHECK_EQUAL(test_1.at(0), 0ul);
            TEST_CHECK_EQUAL(test_1.at(1), 2ul);
            TEST_CHECK_EQUAL(test_1.at(2), 3ul);
            TEST_CHECK_EQUAL(test_1.at(3), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_2(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_edge, 1));
            TEST_CHECK_EQUAL(test_2.size(), 4ul);
            TEST_CHECK_EQUAL(test_2.at(0), 1ul);
            TEST_CHECK_EQUAL(test_2.at(1), 3ul);
            TEST_CHECK_EQUAL(test_2.at(2), 4ul);
            TEST_CHECK_EQUAL(test_2.at(3), 6ul);

            //testing face-face access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_3(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_face, 0));
            TEST_CHECK_EQUAL(test_3.size(), 2ul);
            TEST_CHECK_EQUAL(test_3.at(0), 0ul);
            TEST_CHECK_EQUAL(test_3.at(1), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_4(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_face, 1));
            TEST_CHECK_EQUAL(test_4.size(), 2ul);
            TEST_CHECK_EQUAL(test_4.at(0), 1ul);
            TEST_CHECK_EQUAL(test_4.at(1), 0ul);

            //testing face-vertex access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_5(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_vertex, 0));
            TEST_CHECK_EQUAL(test_5.size(), 4ul);
            TEST_CHECK_EQUAL(test_5.at(0), 0ul);
            TEST_CHECK_EQUAL(test_5.at(1), 1ul);
            TEST_CHECK_EQUAL(test_5.at(2), 3ul);
            TEST_CHECK_EQUAL(test_5.at(3), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_6(m3.get_adjacent_polytopes(fem::pl_face, fem::pl_vertex, 1));
            TEST_CHECK_EQUAL(test_6.size(), 4ul);
            TEST_CHECK_EQUAL(test_6.at(0), 1ul);
            TEST_CHECK_EQUAL(test_6.at(1), 2ul);
            TEST_CHECK_EQUAL(test_6.at(2), 4ul);
            TEST_CHECK_EQUAL(test_6.at(3), 5ul);

            //testing edge-vertex
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_7(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 0));
            TEST_CHECK_EQUAL(test_7.size(), 2ul);
            TEST_CHECK_EQUAL(test_7.at(0), 0ul);
            TEST_CHECK_EQUAL(test_7.at(1), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_8(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 1));
            TEST_CHECK_EQUAL(test_8.size(), 2ul);
            TEST_CHECK_EQUAL(test_8.at(0), 1ul);
            TEST_CHECK_EQUAL(test_8.at(1), 2ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_9(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 2));
            TEST_CHECK_EQUAL(test_9.size(), 2ul);
            TEST_CHECK_EQUAL(test_9.at(0), 0ul);
            TEST_CHECK_EQUAL(test_9.at(1), 3ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_10(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 3));
            TEST_CHECK_EQUAL(test_10.size(), 2ul);
            TEST_CHECK_EQUAL(test_10.at(0), 1ul);
            TEST_CHECK_EQUAL(test_10.at(1), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_11(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 4));
            TEST_CHECK_EQUAL(test_11.size(), 2ul);
            TEST_CHECK_EQUAL(test_11.at(0), 2ul);
            TEST_CHECK_EQUAL(test_11.at(1), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_12(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 5));
            TEST_CHECK_EQUAL(test_12.size(), 2ul);
            TEST_CHECK_EQUAL(test_12.at(0), 3ul);
            TEST_CHECK_EQUAL(test_12.at(1), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_13(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_vertex, 6));
            TEST_CHECK_EQUAL(test_13.size(), 2ul);
            TEST_CHECK_EQUAL(test_13.at(0), 4ul);
            TEST_CHECK_EQUAL(test_13.at(1), 5ul);

            //testing edge-edge
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_14(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 0));
            TEST_CHECK_EQUAL(test_14.size(), 4ul);
            TEST_CHECK_EQUAL(test_14.at(0), 0ul);
            TEST_CHECK_EQUAL(test_14.at(1), 2ul);
            TEST_CHECK_EQUAL(test_14.at(2), 1ul);
            TEST_CHECK_EQUAL(test_14.at(3), 3ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_15(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 1));
            TEST_CHECK_EQUAL(test_15.size(), 4ul);
            TEST_CHECK_EQUAL(test_15.at(0), 0ul);
            TEST_CHECK_EQUAL(test_15.at(1), 1ul);
            TEST_CHECK_EQUAL(test_15.at(2), 3ul);
            TEST_CHECK_EQUAL(test_15.at(3), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_16(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 2));
            TEST_CHECK_EQUAL(test_16.size(), 3ul);
            TEST_CHECK_EQUAL(test_16.at(0), 0ul);
            TEST_CHECK_EQUAL(test_16.at(1), 2ul);
            TEST_CHECK_EQUAL(test_16.at(2), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_17(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 3));
            TEST_CHECK_EQUAL(test_17.size(), 5ul);
            TEST_CHECK_EQUAL(test_17.at(0), 0ul);
            TEST_CHECK_EQUAL(test_17.at(1), 1ul);
            TEST_CHECK_EQUAL(test_17.at(2), 3ul);
            TEST_CHECK_EQUAL(test_17.at(3), 5ul);
            TEST_CHECK_EQUAL(test_17.at(4), 6ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_18(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 4));
            TEST_CHECK_EQUAL(test_18.size(), 3ul);
            TEST_CHECK_EQUAL(test_18.at(0), 1ul);
            TEST_CHECK_EQUAL(test_18.at(1), 4ul);
            TEST_CHECK_EQUAL(test_18.at(2), 6ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_19(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 5));
            TEST_CHECK_EQUAL(test_19.size(), 4ul);
            TEST_CHECK_EQUAL(test_19.at(0), 2ul);
            TEST_CHECK_EQUAL(test_19.at(1), 5ul);
            TEST_CHECK_EQUAL(test_19.at(2), 3ul);
            TEST_CHECK_EQUAL(test_19.at(3), 6ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_20(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_edge, 6));
            TEST_CHECK_EQUAL(test_20.size(), 4ul);
            TEST_CHECK_EQUAL(test_20.at(0), 3ul);
            TEST_CHECK_EQUAL(test_20.at(1), 5ul);
            TEST_CHECK_EQUAL(test_20.at(2), 6ul);
            TEST_CHECK_EQUAL(test_20.at(3), 4ul);

            //testing edge-face access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_21(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 0));
            TEST_CHECK_EQUAL(test_21.size(), 1ul);
            TEST_CHECK_EQUAL(test_21.at(0), 0ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_22(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 1));
            TEST_CHECK_EQUAL(test_22.size(), 1ul);
            TEST_CHECK_EQUAL(test_22.at(0), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_23(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 2));
            TEST_CHECK_EQUAL(test_23.size(), 1ul);
            TEST_CHECK_EQUAL(test_23.at(0), 0ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_24(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 3));
            TEST_CHECK_EQUAL(test_24.size(), 2ul);
            TEST_CHECK_EQUAL(test_24.at(0), 0ul);
            TEST_CHECK_EQUAL(test_24.at(1), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_25(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 4));
            TEST_CHECK_EQUAL(test_25.size(), 1ul);
            TEST_CHECK_EQUAL(test_25.at(0), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_26(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 5));
            TEST_CHECK_EQUAL(test_26.size(), 1ul);
            TEST_CHECK_EQUAL(test_26.at(0), 0ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_27(m3.get_adjacent_polytopes(fem::pl_edge, fem::pl_face, 6));
            TEST_CHECK_EQUAL(test_27.size(), 1ul);
            TEST_CHECK_EQUAL(test_27.at(0), 1ul);

            //testing vertex-vertex access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_28(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 0));
            TEST_CHECK_EQUAL(test_28.size(), 3ul);
            TEST_CHECK_EQUAL(test_28.at(0), 0ul);
            TEST_CHECK_EQUAL(test_28.at(1), 1ul);
            TEST_CHECK_EQUAL(test_28.at(2), 3ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_29(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 1));
            TEST_CHECK_EQUAL(test_29.size(), 4ul);
            TEST_CHECK_EQUAL(test_29.at(0), 0ul);
            TEST_CHECK_EQUAL(test_29.at(1), 1ul);
            TEST_CHECK_EQUAL(test_29.at(2), 2ul);
            TEST_CHECK_EQUAL(test_29.at(3), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_30(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 2));
            TEST_CHECK_EQUAL(test_30.size(), 3ul);
            TEST_CHECK_EQUAL(test_30.at(0), 1ul);
            TEST_CHECK_EQUAL(test_30.at(1), 2ul);
            TEST_CHECK_EQUAL(test_30.at(2), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_31(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 3));
            TEST_CHECK_EQUAL(test_31.size(), 3ul);
            TEST_CHECK_EQUAL(test_31.at(0), 0ul);
            TEST_CHECK_EQUAL(test_31.at(1), 3ul);
            TEST_CHECK_EQUAL(test_31.at(2), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_32(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 4));
            TEST_CHECK_EQUAL(test_32.size(), 4ul);
            TEST_CHECK_EQUAL(test_32.at(0), 1ul);
            TEST_CHECK_EQUAL(test_32.at(1), 4ul);
            TEST_CHECK_EQUAL(test_32.at(2), 3ul);
            TEST_CHECK_EQUAL(test_32.at(3), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_33(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_vertex, 5));
            TEST_CHECK_EQUAL(test_33.size(), 3ul);
            TEST_CHECK_EQUAL(test_33.at(0), 2ul);
            TEST_CHECK_EQUAL(test_33.at(1), 5ul);
            TEST_CHECK_EQUAL(test_33.at(2), 4ul);

            //testing vertex-edge access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_34(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 0));
            TEST_CHECK_EQUAL(test_34.size(), 2ul);
            TEST_CHECK_EQUAL(test_34.at(0), 0ul);
            TEST_CHECK_EQUAL(test_34.at(1), 2ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_35(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 1));
            TEST_CHECK_EQUAL(test_35.size(), 3ul);
            TEST_CHECK_EQUAL(test_35.at(0), 0ul);
            TEST_CHECK_EQUAL(test_35.at(1), 1ul);
            TEST_CHECK_EQUAL(test_35.at(2), 3ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_36(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 2));
            TEST_CHECK_EQUAL(test_36.size(), 2ul);
            TEST_CHECK_EQUAL(test_36.at(0), 1ul);
            TEST_CHECK_EQUAL(test_36.at(1), 4ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_37(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 3));
            TEST_CHECK_EQUAL(test_37.size(), 2ul);
            TEST_CHECK_EQUAL(test_37.at(0), 2ul);
            TEST_CHECK_EQUAL(test_37.at(1), 5ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_38(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 4));
            TEST_CHECK_EQUAL(test_38.size(), 3ul);
            TEST_CHECK_EQUAL(test_38.at(0), 3ul);
            TEST_CHECK_EQUAL(test_38.at(1), 5ul);
            TEST_CHECK_EQUAL(test_38.at(2), 6ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_39(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_edge, 5));
            TEST_CHECK_EQUAL(test_39.size(), 2ul);
            TEST_CHECK_EQUAL(test_39.at(0), 4ul);
            TEST_CHECK_EQUAL(test_39.at(1), 6ul);

            //testing vertex-face access
            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_40(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 0));
            TEST_CHECK_EQUAL(test_40.size(), 1ul);
            TEST_CHECK_EQUAL(test_40.at(0), 0ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_41(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 1));
            TEST_CHECK_EQUAL(test_41.size(), 2ul);
            TEST_CHECK_EQUAL(test_41.at(0), 0ul);
            TEST_CHECK_EQUAL(test_41.at(1), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_42(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 2));
            TEST_CHECK_EQUAL(test_42.size(), 1ul);
            TEST_CHECK_EQUAL(test_42.at(0), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_43(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 3));
            TEST_CHECK_EQUAL(test_43.size(), 1ul);
            TEST_CHECK_EQUAL(test_43.at(0), 0ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_44(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 4));
            TEST_CHECK_EQUAL(test_44.size(), 2ul);
            TEST_CHECK_EQUAL(test_44.at(0), 0ul);
            TEST_CHECK_EQUAL(test_44.at(1), 1ul);

            typename fem::Topology<IndexType_, OT_, IT_>::storage_type_ test_45(m3.get_adjacent_polytopes(fem::pl_vertex, fem::pl_face, 5));
            TEST_CHECK_EQUAL(test_45.size(), 1ul);
            TEST_CHECK_EQUAL(test_45.at(0), 1ul);
        }
};
MeshTest<tags::CPU, unsigned long, std::vector, std::vector<unsigned long> > topology_test_cpu_v_v("std::vector, std::vector");
MeshTest<tags::CPU, unsigned long, std::deque, std::vector<unsigned long> > topology_test_cpu_d_v("std::deque, std::vector");
MeshTest<tags::CPU, unsigned long, std::vector, std::deque<unsigned long> > topology_test_cpu_v_d("std::vector, std::deque");
MeshTest<tags::CPU, unsigned long, std::deque, std::deque<unsigned long> > topology_test_cpu_d_d("std::deque, std::deque");

//MeshTest<tags::CPU, unsigned long, std::vector, fem::DenseDataWrapper<15, unsigned long> > topology_test_cpu_v_ddw("std::vector, DV");
