/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef FEM_GUARD_DENSE_DATA_WRAPPER_HH
#define FEM_GUARD_DENSE_DATA_WRAPPER_HH 1

#include <honei/la/dense_vector.hh>

namespace honei
{
    namespace fem
    {
        class DenseVector<unsigned long>;
        template<unsigned long _i, typename DT_, template<typename> class ContType_ = DenseVector>
        class DenseDataWrapper
        {
            public:
                DenseDataWrapper() :
                    _size(_i),
                    _num_non_zeros(0),
                    _data(new ContType_<DT_>(_i))
                {
                }

                ~DenseDataWrapper()
                {
                    delete _data;
                }

                unsigned long size()
                {
                    return _num_non_zeros;
                }

                unsigned long capacity()
                {
                    return _size - _num_non_zeros;
                }

                void push_back(DT_ value)
                {
                    //todo capacity check
                    (*_data)[_num_non_zeros] = value;
                    ++_num_non_zeros;
                }

                DT_ & at(unsigned long i)
                {
                    //todo in non-zero range check
                    return (*_data)[i];
                }

                DT_ & operator[](unsigned long i)
                {
                    //todo in non-zero range check
                    return (*_data)[i];
                }

            private:
                unsigned long _size;
                unsigned long _num_non_zeros;
                ContType_<DT_>* _data;
        };
    }
}
#endif
