/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_MESH_HH
#define FEM_GUARD_MESH_HH 1

#include <honei/fem/topology.hh>

namespace honei
{
    namespace fem
    {
        enum RequiredNumTopologies
        {
            rnt_1D = 2,
            rnt_2D = 3,
            rnt_3D = 4
        };

        enum PolytopeLevels
        {
            pl_vertex = 0,
            pl_edge,
            pl_face,
            pl_polyhedron
        };

        template<unsigned _i = rnt_2D, typename TopologyType_ = Topology<> >
        class Mesh
        {
            public:
                Mesh() :
                    _num_inter_topologies(_i),
                    _inter_topology_relations(new unsigned[_i]),
                    _inter_topologies(new TopologyType_[_i]),
                    _element_topology(new TopologyType_)
                {
                    const unsigned i_end(_i);
                    for(unsigned i(0) ; i < i_end ; ++i)
                    {
                        _inter_topology_relations[i] = i == 0 ? _i - 1 : i - 1;
                    }

                }

                ~Mesh()
                {
                    delete[] _inter_topologies;
                    delete[] _inter_topology_relations;
                    delete _element_topology;
                }

                unsigned get_inter_topology_relation(unsigned i)
                {
                    return _inter_topology_relations[i];
                }

                template<typename NLT_>
                void add_element(const NLT_ & neighbours)
                {
                    _element_topology->push_back(neighbours);
                }

                void add_element()
                {
                    _element_topology->push_back();
                }

                template<typename NLT_>
                void add_polytope(unsigned level, const NLT_ & neighbours)
                {
                    if(level == _i)
                        add_element(neighbours);
                    else if(level < _i)
                        _inter_topologies[level].push_back(neighbours);
                    //todo else exception
                }

                void add_polytope(unsigned level)
                {
                    if(level == _i)
                        add_element();
                    else if(level < _i)
                        _inter_topologies[level].push_back();
                }

            private:
                unsigned _num_inter_topologies;
                unsigned * _inter_topology_relations;

                TopologyType_* _inter_topologies;
                TopologyType_* _element_topology;

        };
    }
}
#endif
