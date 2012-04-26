/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_MESH_HH
#define FEM_GUARD_MESH_HH 1

#include <honei/fem/topology.hh>

namespace honei
{
    namespace fem
    {
        template<typename TopologyType_ = Topology<> >
        class Mesh
        {
            public:
                Mesh(unsigned long i) :
                    _num_inter_topologies(i),
                    _inter_topologies(new TopologyType_[i]),
                    _element_topology(new TopologyType_)
                {
                }

                ~Mesh()
                {
                    delete[] _inter_topologies;
                    delete _element_topology;
                }

            private:
                unsigned long _num_inter_topologies;
                TopologyType_* _inter_topologies;
                TopologyType_* _element_topology;
        };
    }
}
#endif
