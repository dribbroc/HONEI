/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_HALO_HH
#define FEM_GUARD_HALO_HH 1

#include<vector>

namespace honei
{
    namespace fem
    {
        template<unsigned delta, typename MeshType_, template<typename, typename> class StorageType_ = std::vector, typename IndexType_ = unsigned long>
        class Halo
        {
            public:
                typedef IndexType_ index_type_;
                typedef MeshType_ mesh_type_;

                Halo(MeshType_ & left, MeshType_ & right) :
                    _halo_elements(StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _halo_element_counterparts(StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _left(left),
                    _right(right),
                    _overlap(delta)
                {
                }

                ~Halo()
                {
                }

                ///Add correspondence of i
                void add_halo_element_pair(IndexType_ i, IndexType_ j)
                {
                    _halo_elements.push_back(i);
                    _halo_element_counterparts.push_back(j);
                }

                IndexType_ get_element_counterpart(IndexType_ index)
                {
                    return _halo_element_counterparts.at(index);
                }

                IndexType_ get_element(IndexType_ index)
                {
                    return _halo_elements.at(index);
                }

                IndexType_ size()
                {
                    return _halo_elements.size();
                }

                MeshType_ & get_left_mesh()
                {
                    return _left;
                }

                MeshType_ & get_right_mesh()
                {
                    return _right;
                }

                unsigned get_overlap()
                {
                    return _overlap;
                }

            private:
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_elements;
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_element_counterparts;

                MeshType_ & _left;
                MeshType_ & _right;

                unsigned _overlap;
        };

    }

}


#endif
