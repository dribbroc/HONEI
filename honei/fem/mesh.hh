/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_MESH_HH
#define FEM_GUARD_MESH_HH 1

#include <honei/fem/topology.hh>
#include <honei/fem/mesh_error.hh>
#include <iostream>
#include <cmath>

namespace honei
{
    namespace fem
    {
        template<typename MeshType_, typename ExplicitType_>
        class MeshAttributeRegistration
        {
        };

        template<typename MeshType_>
        class MeshAttributeRegistration<MeshType_, typename MeshType_::attr_type_1_>
        {
            public:
                static unsigned execute(MeshType_ & mesh, const unsigned polytope_level)
                {
                    typename MeshType_::attribute_storage_type_1_ storage;
                    mesh._attributes_of_type_1->push_back(storage);
                    mesh._attribute_polytopelevel_relations_1->push_back(polytope_level);

                    ++mesh._num_attributes_of_type_1;
                    return mesh._num_attributes_of_type_1 - 1;
                }
        };

        template<typename MeshType_>
        class MeshAttributeRegistration<MeshType_, typename MeshType_::attr_type_2_>
        {
            public:
                static unsigned execute(MeshType_ & mesh, const unsigned polytope_level)
                {
                    typename MeshType_::attribute_storage_type_2_ storage;
                    mesh._attributes_of_type_2->push_back(storage);
                    mesh._attribute_polytopelevel_relations_2->push_back(polytope_level);

                    ++mesh._num_attributes_of_type_2;
                    return mesh._num_attributes_of_type_2 - 1;
                }
        };

        template<typename MeshType_>
        class MeshAttributeRegistration<MeshType_, typename MeshType_::attr_type_3_>
        {
            public:
                static unsigned execute(MeshType_ & mesh, const unsigned polytope_level)
                {
                    typename MeshType_::attribute_storage_type_2_ storage;
                    mesh._attributes_of_type_2->push_back(storage);
                    mesh._attribute_polytopelevel_relations_2->push_back(polytope_level);

                    ++mesh._num_attributes_of_type_2;
                    return mesh._num_attributes_of_type_2 - 1;
                }
        };

        enum RequiredNumTopologies
        {
            rnt_1D = 2,
            rnt_2D = 4,
            rnt_3D = 6
        };

        enum PolytopeLevels
        {
            pl_vertex = 0,
            pl_edge,
            pl_face,
            pl_polyhedron
        };

        enum InternalPolytopeIndices
        {
            ipi_vertex_edge = 0,
            ipi_edge_vertex,
            ipi_edge_face,
            ipi_face_edge,
            ipi_face_polyhedron,
            ipi_polyhedron_face
        };

        template<RequiredNumTopologies _i = rnt_2D,
                 typename TopologyType_ = Topology<>,
                 template<typename, typename> class AttributeStorageType_ = std::vector,
                 template<typename, typename> class OuterAttributeStorageType_ = std::vector,
                 typename AttributeType1_ = double,
                 typename AttributeType2_ = unsigned long,
                 typename AttributeType3_ = float>
        class Mesh
        {
            public:
                friend class MeshAttributeRegistration< Mesh<_i, TopologyType_, AttributeStorageType_, OuterAttributeStorageType_, AttributeType1_, AttributeType2_>, AttributeType1_>;
                friend class MeshAttributeRegistration< Mesh<_i, TopologyType_, AttributeStorageType_, OuterAttributeStorageType_, AttributeType1_, AttributeType2_>, AttributeType2_>;
                friend class MeshAttributeRegistration< Mesh<_i, TopologyType_, AttributeStorageType_, OuterAttributeStorageType_, AttributeType1_, AttributeType3_>, AttributeType3_>;

                typedef AttributeType1_ attr_type_1_;
                typedef AttributeType2_ attr_type_2_;

                typedef OuterAttributeStorageType_<
                    AttributeStorageType_<
                        AttributeType1_,
                        std::allocator<AttributeType1_> >,
                    std::allocator<
                            AttributeStorageType_<
                                AttributeType1_,
                                std::allocator<AttributeType1_>
                > > > outer_attribute_storage_type_1_;

                typedef OuterAttributeStorageType_<
                    AttributeStorageType_<
                        AttributeType2_,
                        std::allocator<AttributeType2_> >,
                    std::allocator<
                            AttributeStorageType_<
                                AttributeType2_,
                                std::allocator<AttributeType2_>
                > > > outer_attribute_storage_type_2_;

                typedef OuterAttributeStorageType_<
                    AttributeStorageType_<
                        AttributeType2_,
                        std::allocator<AttributeType2_> >,
                    std::allocator<
                            AttributeStorageType_<
                                AttributeType2_,
                                std::allocator<AttributeType2_>
                > > > outer_attribute_storage_type_3_;

                typedef AttributeStorageType_<
                            AttributeType1_,
                            std::allocator<AttributeType1_> > attribute_storage_type_1_;

                typedef AttributeStorageType_<
                            AttributeType2_,
                            std::allocator<AttributeType2_> > attribute_storage_type_2_;

                typedef AttributeStorageType_<
                            AttributeType2_,
                            std::allocator<AttributeType3_> > attribute_storage_type_3_;

                Mesh() :
                    _num_inter_topologies(_i),
                    _num_levels((unsigned)(_i/2u) + 1u),
                    _topologies(new TopologyType_[_i]),
                    _num_attributes_of_type_1(0),
                    _num_attributes_of_type_2(0),
                    _num_attributes_of_type_3(0),
                    _attribute_polytopelevel_relations_1(new typename TopologyType_::storage_type_),
                    _attribute_polytopelevel_relations_2(new typename TopologyType_::storage_type_),
                    _attribute_polytopelevel_relations_3(new typename TopologyType_::storage_type_),
                    _attributes_of_type_1(new outer_attribute_storage_type_1_),
                    _attributes_of_type_2(new outer_attribute_storage_type_2_),
                    _attributes_of_type_3(new outer_attribute_storage_type_3_)
                {
                }

                ~Mesh()
                {
                    delete[] _topologies;
                }

                void add_polytope(const unsigned level)
                {
                    CONTEXT("When adding polytope");
#ifdef FEM_MESH_DEBUG
                    if(level > _num_levels - 1)
                        throw MeshInternalIndexOutOfBounds(level, _num_levels - 1);
#endif
                    if(_upward_index(level) != -1)
                        _topologies[_upward_index(level)].push_back();

                    if(_downward_index(level) != -1)
                        _topologies[_downward_index(level)].push_back();
                }

                template<typename IndexT_>
                void add_adjacency(const unsigned from_level,
                                   const unsigned to_level,
                                   const IndexT_ polytope_index,
                                   const IndexT_ value)
                {
                    CONTEXT("When adding adjacency");
#ifdef FEM_MESH_DEBUG
                    if(_level_difference(from_level, to_level) > 1)
                        throw MeshError("Level difference too large.");

                    if(from_level > _num_levels - 1)
                        throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);

                    if(to_level > _num_levels - 1)
                        throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);

                    if(from_level == to_level)
                        throw MeshError("Explicit neighbourhood of polytopes on same level currently unsupported. Use implict storage instead.");//TODO
#endif

                    if(from_level < to_level)
                    {
#ifdef FEM_MESH_DEBUG
                        if(_upward_index(from_level) == -1)
                            throw MeshError("Invalid impairment of levels.");
#endif
                        _topologies[_upward_index(from_level)].at(polytope_index).push_back(value);

                        if(_downward_index(to_level) != -1)
                            _topologies[_downward_index(to_level)].at(value).push_back(polytope_index);
                    }
                    else if(from_level > to_level)
                    {
#ifdef FEM_MESH_DEBUG
                        if(_downward_index(from_level) == -1)
                            throw MeshError("Invalid impairment of levels.");
#endif
                        _topologies[_downward_index(from_level)].at(polytope_index).push_back(value);

                        if(_upward_index(to_level) != -1)
                            _topologies[_upward_index(to_level)].at(value).push_back(polytope_index);
                    }
                }

                template<typename IndexType_>
                typename TopologyType_::storage_type_ get_adjacent_polytopes(const unsigned from_level, const unsigned to_level, IndexType_ i)
                {
                    CONTEXT("When calculating adjacent polytopes");

#ifdef FEM_MESH_DEBUG
                    if(from_level > _num_levels - 1)
                        throw MeshInternalIndexOutOfBounds(from_level, _num_levels - 1);
                    if(to_level > _num_levels - 1)
                        throw MeshInternalIndexOutOfBounds(to_level, _num_levels - 1);
#endif

                    const int sweepdir(_sweep_direction(from_level, to_level));
                    int current_sweepdir(sweepdir == 0 ? -1 : sweepdir);
                    current_sweepdir = from_level == 0 && sweepdir == 0 ? 1 : current_sweepdir;
                    const unsigned level_diff(_level_difference(from_level, to_level));
                    const unsigned path_length(level_diff + 1);
                    unsigned current_level(from_level);
                    typename TopologyType_::compound_storage_type_ search_data;


                    if(level_diff == 0u) //we can directly give the neighbours
                    {
                        if(sweepdir == -1) //sweep down
                            return _topologies[_downward_index(from_level)].at(i);
                        else if(sweepdir == 1) //sweep up
                            return _topologies[_upward_index(from_level)].at(i);
#ifdef FEM_MESH_DEBUG
                        else
                            throw MeshError("Sweep direction / sweep length mismatch.");
#endif
                    }
                    else
                    {

                        for(unsigned j(0) ; j < path_length ; ++j)
                        {
                            //create current search datastructure
                            typename TopologyType_::storage_type_ temp;
                            search_data.push_back(temp);

                            //fill current datastructure:
                            if(j == 0)
                            {
                                search_data.at(0) = current_sweepdir == -1 ? _topologies[_downward_index(current_level)].at(i)
                                                                           : _topologies[_upward_index(current_level)].at(i);

                                current_level = current_sweepdir == 1 ? current_level + 1
                                                                      : (current_sweepdir == 0 ? current_level + 1
                                                                                               : current_level - 1);


                                current_sweepdir = sweepdir == 0 && current_sweepdir == -1 ? 1
                                                                                           : (sweepdir == 0 && current_sweepdir == 1 ? -1
                                                                                                                                     : current_sweepdir);
                            }
                            else
                            {
                                //for all entries in previous polytope list get all polytope lists and store sequentially
                                for(IndexType_ k(0) ; k < (IndexType_)search_data.at(j - 1).size(); ++k)
                                {
                                    IndexType_ l_upper(
                                            current_sweepdir == -1 ? (IndexType_)_topologies[_downward_index(current_level)].at(search_data.at(j - 1).at(k)).size()
                                                                   : (IndexType_)_topologies[_upward_index(current_level)].at(search_data.at(j - 1).at(k)).size()
                                                      );

                                    for(IndexType_ l(0) ; l < l_upper ; ++l)
                                    {
                                        //TODO optimise search
                                        IndexType_ to_insert(
                                                current_sweepdir == -1 ? _topologies[_downward_index(current_level)].at(search_data.at(j - 1).at(k)).at(l)
                                                                       : _topologies[_upward_index(current_level)].at(search_data.at(j - 1).at(k)).at(l)
                                                            );

                                        bool insert(true);
                                        for(IndexType_ search_index(0) ; search_index < (IndexType_)search_data.at(j).size() ; ++search_index)
                                        {
                                            if((IndexType_)search_data.at(j).at(search_index) == (IndexType_)to_insert)
                                            {
                                                insert = false;
                                                break;
                                            }
                                        }
                                        if(insert)
                                        {
                                            search_data.at(j).push_back(to_insert);
                                        }
                                    }
                                }
                                current_level = current_sweepdir == 1 ? current_level + 1
                                                                      : current_level - 1;
                            }
                        }
                        return search_data.at(search_data.size() - 1);
                    }
                    return search_data.at(search_data.size() - 1);
                }

                const typename TopologyType_::index_type_ get_num_levels()
                {
                    return _num_levels;
                }

                const int get_downward_index(const unsigned pl)
                {
                    return _downward_index(pl);
                }

                const int get_upward_index(const unsigned pl)
                {
                    return _upward_index(pl);
                }

                template<typename AT_>
                void set_attribute_value(const unsigned attribute_index, const typename TopologyType_::index_type_ index, const AT_ value)
                {
                    if(typeid(AT_) == typeid(AttributeType1_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_1 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_1)
                        {
                            _attributes_of_type_1->at(attribute_index).at(index) = value;
                        }
                    }
                    else if(typeid(AT_) == typeid(AttributeType2_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_2 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_2)
                        {
                            _attributes_of_type_2->at(attribute_index).at(index) = value;
                        }
                    }
                    else if(typeid(AT_) == typeid(AttributeType3_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_3 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_3)
                        {
                            _attributes_of_type_3->at(attribute_index).at(index) = value;
                        }
                    }
                    //todo catch index out of bounds
                }

                template<typename AT_>
                void add_attribute_value(const unsigned attribute_index, const AT_ value)
                {
                    if(typeid(AT_) == typeid(AttributeType1_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_1 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_1)
                        {
                            _attributes_of_type_1->at(attribute_index).push_back(value);
                        }
                    }
                    else if(typeid(AT_) == typeid(AttributeType2_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_2 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_2)
                        {
                            _attributes_of_type_2->at(attribute_index).push_back(value);
                        }
                    }
                    else if(typeid(AT_) == typeid(AttributeType3_))
                    {
#ifdef FEM_MESH_DEBUG
                        if(_num_attributes_of_type_3 == 0)
                            throw MeshInternalIndexOutOfBounds(attribute_index, 0);
#endif
                        if(attribute_index < _num_attributes_of_type_3)
                        {
                            _attributes_of_type_3->at(attribute_index).push_back(value);
                        }
                    }
                }

                template<typename IndexType_>
                typename TopologyType_::storage_type_ get_primary_comm_neighbours(IndexType_ i)
                {
                    typename TopologyType_::storage_type_ pre_result(this->get_adjacent_polytopes(_num_levels - 1, _num_levels - 1, i));
                    typename TopologyType_::storage_type_ result;

                    //remove self
                    for(IndexType_ j(0) ; j < (IndexType_)pre_result.size() ; ++j)
                        if((IndexType_)pre_result.at(j) != i)
                            result.push_back(pre_result.at(j));

                    return result;
                }

                template<typename IndexType_>
                typename TopologyType_::storage_type_ get_all_comm_neighbours(IndexType_ i)
                {
                    const int sweepdir(-1);
                    int current_sweepdir(sweepdir);
                    const unsigned path_length(_num_inter_topologies);
                    unsigned current_level(_num_levels - 1);
                    typename TopologyType_::compound_storage_type_ search_data;

                    for(unsigned j(0) ; j < path_length ; ++j)
                    {
                        //create current search datastructure
                        typename TopologyType_::storage_type_ temp;
                        search_data.push_back(temp);

                        //fill current datastructure:
                        if(j == 0)
                        {
                            search_data.at(0) = _topologies[_downward_index(current_level)].at(i);

                            --current_level;

                            current_sweepdir = current_level == 0 ? 1 : current_sweepdir;
                        }
                        else
                        {
                            //for all entries in previous polytope list get all polytope lists and store sequentially
                            for(IndexType_ k(0) ; k < (IndexType_)search_data.at(j - 1).size(); ++k)
                            {
                                IndexType_ l_upper(
                                        current_sweepdir == -1 ? (IndexType_)_topologies[_downward_index(current_level)].at(search_data.at(j - 1).at(k)).size()
                                                               : (IndexType_)_topologies[_upward_index(current_level)].at(search_data.at(j - 1).at(k)).size()
                                        );

                                for(IndexType_ l(0) ; l < l_upper ; ++l)
                                {
                                    //TODO optimise search
                                    IndexType_ to_insert(
                                            current_sweepdir == -1 ? _topologies[_downward_index(current_level)].at(search_data.at(j - 1).at(k)).at(l)
                                                                   : _topologies[_upward_index(current_level)].at(search_data.at(j - 1).at(k)).at(l)
                                            );

                                    bool insert(true);
                                    for(IndexType_ search_index(0) ; search_index < (IndexType_)search_data.at(j).size() ; ++search_index)
                                    {
                                        if((IndexType_)search_data.at(j).at(search_index) == (IndexType_)to_insert)
                                        {
                                            insert = false;
                                            break;
                                        }
                                    }
                                    if(insert)
                                    {
                                        search_data.at(j).push_back(to_insert);
                                    }
                                }
                            }
                            current_level = current_sweepdir == 1 ? current_level + 1
                                                                  : current_level - 1;

                            current_sweepdir = current_level == 0 ? 1 : current_sweepdir;
                        }
                    }

                    typename TopologyType_::storage_type_ result;

                    //remove self
                    for(IndexType_ j(0) ; j < (IndexType_)search_data.at(search_data.size() - 1).size() ; ++j)
                        if((IndexType_)search_data.at(search_data.size() - 1).at(j) != i)
                            result.push_back(search_data.at(search_data.size() - 1).at(j));

                    return result;
                }

                unsigned get_num_attributes_of_type_1()
                {
                    return _num_attributes_of_type_1;
                }

                unsigned get_num_attributes_of_type_2()
                {
                    return _num_attributes_of_type_2;
                }

                unsigned get_num_attributes_of_type_3()
                {
                    return _num_attributes_of_type_3;
                }

                typename TopologyType_::storage_type_ & get_attribute_polytopelevel_relations_1()
                {
                    return *_attribute_polytopelevel_relations_1;
                }

                typename TopologyType_::storage_type_ & get_attribute_polytopelevel_relations_2()
                {
                    return *_attribute_polytopelevel_relations_2;
                }

                typename TopologyType_::storage_type_ & get_attribute_polytopelevel_relations_3()
                {
                    return *_attribute_polytopelevel_relations_3;
                }

                outer_attribute_storage_type_1_ & get_attributes_of_type_1()
                {
                    return *_attributes_of_type_1;
                }

                outer_attribute_storage_type_2_ & get_attributes_of_type_2()
                {
                    return *_attributes_of_type_2;
                }

                outer_attribute_storage_type_3_ & get_attributes_of_type_3()
                {
                    return *_attributes_of_type_3;
                }

            private:
                const unsigned _num_inter_topologies;
                const unsigned _num_levels;
                TopologyType_* _topologies;

                unsigned _num_attributes_of_type_1;
                unsigned _num_attributes_of_type_2;
                unsigned _num_attributes_of_type_3;

                typename TopologyType_::storage_type_ * _attribute_polytopelevel_relations_1;
                typename TopologyType_::storage_type_ * _attribute_polytopelevel_relations_2;
                typename TopologyType_::storage_type_ * _attribute_polytopelevel_relations_3;

                outer_attribute_storage_type_1_* _attributes_of_type_1;
                outer_attribute_storage_type_2_* _attributes_of_type_2;
                outer_attribute_storage_type_3_* _attributes_of_type_3;

                inline const unsigned _level_difference(const unsigned from, const unsigned to)
                {
                    return from == to ? 1u : ( from > to ? (unsigned)std::abs((double)(from - to)) - 1u : (unsigned)std::abs((double)(to - from)) - 1u);
                }

                inline const int _downward_index(const unsigned pl)
                {
                    return (pl == 0u || pl >= _num_levels) ? - 1 : (pl == 3u ? 5 : (unsigned)pow(2, pl) - 1);
                }

                inline const int _upward_index(const unsigned pl)
                {
                    return pl >= _num_levels - 1? -1 : ( pl > 0 ? (unsigned)pow(2, pl) : 0);
                }

                inline const int _sweep_direction(const unsigned from_level, const unsigned to_level)
                {
                    return from_level == to_level ? 0 : (from_level > to_level ? -1 : 1);
                }
        };
    }
}
#endif
