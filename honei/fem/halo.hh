/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_HALO_HH
#define FEM_GUARD_HALO_HH 1

#include<vector>

namespace honei
{
    namespace fem
    {
        template<unsigned delta = 1, template<typename, typename> class StorageType_ = std::vector, typename IndexType_ = unsigned long>
        class Halo
        {
            public:
                typedef IndexType_ index_type_;

                Halo() :
                    _halo_elements(StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _halo_element_counterparts(StorageType_<IndexType_, std::allocator<IndexType_> >())
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

            private:
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_elements;
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_element_counterparts;
        };

        template<template<typename, typename> class StorageType_, typename IndexType_>
        class Halo<0u, StorageType_, IndexType_>
        {
            public:
                typedef IndexType_ index_type_;

                Halo() :
                    _halo_elements(StorageType_<IndexType_, std::allocator<IndexType_> >()),
                    _halo_element_counterparts(StorageType_<IndexType_, std::allocator<IndexType_> >())
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

            private:
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_elements;
                StorageType_<IndexType_, std::allocator<IndexType_> > _halo_element_counterparts;
        };

        /*template<typename HaloType_, typename MeshType_, typename IndexType_>
        static IndexType_ get_element_counterpart(MeshType_ & source, MeshType_ & target, HaloType_ & halo, IndexType_ index)
        {
        }

        template<typename MeshType_, typename IndexType_>
        static IndexType_  get_element_counterpart(MeshType_ & source, MeshType_ & target, Halo<0, std::vector> & halo, IndexType_ index)
        {
            return 48;
        }*/

        /*template<typename IndexType_>
        static StorageType_ & get_primary_comm_neighbours(MeshType_ & m1, MeshType_ & m2, Halo<0, StorageType_, IndexType_> & halo, IndexType_ index)
        {
                    //we know it must refer to targets max polytope_level
                    unsigned element_level(target.get_num_levels() - 1);
        }

        template<typename IndexType_>
        static StorageType_ & get_all_comm_neighbours(MeshType_ & m1, MeshType_ & m2, Halo<0, StorageType_, IndexType_> & halo, IndexType_ index)
        {
                    //we know it must refer to targets max polytope_level
                    unsigned element_level(target.get_num_levels() - 1);
        }*/
    }

}


#endif
