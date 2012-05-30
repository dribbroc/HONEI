/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef QT_GUARD_VERTEX_HH
#define QT_GUARD_VERTEX_HH 1

template<typename DT_>
struct Vertex
{
    DT_ coord_x, coord_y, coord_z;

    Vertex(DT_ x, DT_ y, DT_ z)
        :
            coord_x(x),
            coord_y(y),
            coord_z(z)
    {
    }

    Vertex()
        :
            coord_x(DT_(0.)),
            coord_y(DT_(0.)),
            coord_z(DT_(0.))
    {
    }
};

#endif
