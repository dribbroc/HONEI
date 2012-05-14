/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_MESH_DATA_HH
#define FEM_GUARD_MESH_DATA_HH 1

#include <honei/fem/mesh.hh>
#include <honei/la/dense_vector.hh>
#include <iostream>
#include <cmath>

namespace honei
{
    namespace fem
    {
        template<
            typename MeshType_,
            template<typename, typename> class OuterStorageType_ = std::vector,
            template<typename>  class VectorType_ = DenseVector,
            typename IndexType_ = unsigned long,
            typename DataType1_ = double,
            typename DataType2_ = unsigned long,
            typename DataType3_ = float>
        class MeshData
        {
            public:
                MeshData(MeshType_ & mesh) :
                    _mesh(mesh),
                    _num_attributes_of_type_1(mesh.get_num_attributes_of_type_1()),
                    _num_attributes_of_type_2(mesh.get_num_attributes_of_type_2()),
                    _num_attributes_of_type_3(mesh.get_num_attributes_of_type_3()),
                    _attributes_of_type_1(OuterStorageType_<VectorType_<DataType1_>, std::allocator<VectorType_<DataType1_> > >()),
                    _attributes_of_type_2(OuterStorageType_<VectorType_<DataType2_>, std::allocator<VectorType_<DataType2_> > >()),
                    _attributes_of_type_3(OuterStorageType_<VectorType_<DataType3_>, std::allocator<VectorType_<DataType3_> > >()),
                    _attribute_polytopelevel_relations_1(VectorType_<IndexType_>(_num_attributes_of_type_1)),
                    _attribute_polytopelevel_relations_2(VectorType_<IndexType_>(_num_attributes_of_type_2)),
                    _attribute_polytopelevel_relations_3(VectorType_<IndexType_>(_num_attributes_of_type_3))
                {
                    for(IndexType_ i(0) ; i < _num_attributes_of_type_1 ; ++i)
                    {
                        VectorType_<DataType1_> attr_i((IndexType_)_mesh.get_attributes_of_type_1().at(i).size());

                        _attribute_polytopelevel_relations_1[i] = (IndexType_)_mesh.get_attribute_polytopelevel_relations_1().at(i);

                        for(unsigned long j(0) ; j < (IndexType_)_mesh.get_attributes_of_type_1().at(i).size() ; ++j)
                        {
                            attr_i[j] = _mesh.get_attributes_of_type_1().at(i).at(j);
                        }

                        _attributes_of_type_1.push_back(attr_i);
                    }

                    for(IndexType_ i(0) ; i < _num_attributes_of_type_2 ; ++i)
                    {
                        VectorType_<DataType2_> attr_i((IndexType_)_mesh.get_attributes_of_type_2().at(i).size());

                        _attribute_polytopelevel_relations_2[i] = (IndexType_)_mesh.get_attribute_polytopelevel_relations_2().at(i);

                        for(unsigned long j(0) ; j < (IndexType_)_mesh.get_attributes_of_type_2().at(i).size() ; ++j)
                        {
                            attr_i[j] = _mesh.get_attributes_of_type_2().at(i).at(j);
                        }

                        _attributes_of_type_2.push_back(attr_i);
                    }

                    for(IndexType_ i(0) ; i < _num_attributes_of_type_3 ; ++i)
                    {
                        VectorType_<DataType3_> attr_i((IndexType_)_mesh.get_attributes_of_type_3().at(i).size());

                        _attribute_polytopelevel_relations_3[i] = (IndexType_)_mesh.get_attribute_polytopelevel_relations_3().at(i);

                        for(unsigned long j(0) ; j < (IndexType_)_mesh.get_attributes_of_type_3().at(i).size() ; ++j)
                        {
                            attr_i[j] = _mesh.get_attributes_of_type_3().at(i).at(j);
                        }

                        _attributes_of_type_3.push_back(attr_i);
                    }
                }

                ~MeshData()
                {
                }

                MeshType_ & get_mesh()
                {
                    return _mesh;
                }

                VectorType_<DataType1_> & get_attribute_of_type_1(IndexType_ index)
                {
                    return _attributes_of_type_1.at(index);
                }

                VectorType_<DataType1_> & get_attribute_of_type_2(IndexType_ index)
                {
                    return _attributes_of_type_2.at(index);
                }

                VectorType_<DataType1_> & get_attribute_of_type_3(IndexType_ index)
                {
                    return _attributes_of_type_3.at(index);
                }

            private:
                MeshType_ & _mesh;

                unsigned _num_attributes_of_type_1;
                unsigned _num_attributes_of_type_2;
                unsigned _num_attributes_of_type_3;

                OuterStorageType_<VectorType_<DataType1_>, std::allocator<VectorType_<DataType1_> > > _attributes_of_type_1;
                OuterStorageType_<VectorType_<DataType2_>, std::allocator<VectorType_<DataType2_> > > _attributes_of_type_2;
                OuterStorageType_<VectorType_<DataType3_>, std::allocator<VectorType_<DataType3_> > > _attributes_of_type_3;

                VectorType_<IndexType_> _attribute_polytopelevel_relations_1;
                VectorType_<IndexType_> _attribute_polytopelevel_relations_2;
                VectorType_<IndexType_> _attribute_polytopelevel_relations_3;
        };
    }
}
#endif
