/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Datatypes.h
 *
 *  Created on: Mar 11, 2010
 *      Author: babrodtk
 */

#ifndef NETCDF_DATATYPES_H_
#define NETCDF_DATATYPES_H_

#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif
#include <iostream>
#include <iomanip>
#include <limits>
#include <honei/util/tr1_boost.hh>

#include <honei/util/instantiation_policy.hh>
#include <honei/util/kpnetcdf_types.hh>

/**
 * Very simple class to simply hold the data we need to represent a heightmap
 */
struct Field {
    float* data;
    unsigned int nx, ny;
    float dx, dy;
    float no_data_value;
    KPInitialConditions::DATA_LOCATION location;

    Field(unsigned int width, unsigned int height) {
        nx = width;
        ny = height;
        data = new float[nx*ny];
        no_data_value = std::numeric_limits<float>::max();
        dx = -1;
        dy = -1;
        location = KPInitialConditions::GRID_CELL_UNKNOWN;
    }
    ~Field() {
        delete [] data;
    }
};

typedef shared_ptr<Field> shared_field_ptr_array3[3];
typedef shared_ptr<Field> shared_field_ptr;

inline std::ostream& operator<<(std::ostream& out, const Field& f) {
    out << "Field [" << f.nx << "x" << f.ny << "]" << "x[" << f.dx << "x" << f.dy << "], no_data_value=" << f.no_data_value << std::endl;

    for (unsigned int i=0; i<f.ny; ++i) {
        for (unsigned int j=0; j<f.nx; ++j) {
            out << std::fixed << std::setprecision(3) << std::setw(6) << f.data[f.nx*i+j] << " ";
        }
        if (i < f.ny-1)
            out << std::endl;
    }
    out << std::setprecision(-1);
    return out;
}

inline std::ostream& operator<<(std::ostream& out, const Field* f) {
    out << *f;
    return out;
}

/**
 * Struct that describes a iteration/timestep, e.g., from a previous simulation.
 */
class TimeStep :
        public honei::InstantiationPolicy<TimeStep, honei::NonCopyable>
{
    public:
        TimeStep(unsigned int nx, unsigned int ny) {
            U[0].reset(new Field(nx, ny));
            U[1].reset(new Field(nx, ny));
            U[2].reset(new Field(nx, ny));
            this->nx = nx;
            this->ny = ny;
        }

        /**
         * Water elevation (w) and discharges (hu, hv) given at grid cell
         * corners and grid cell interfaces, respectively.
         */
        shared_field_ptr_array3 U;

        /**
         * Number of nodes for water elevation and discharges
         * The simulator then creates a simulation domain of nx-1 x ny-1 cells.
         * The simulator further uses two ghost cells, making the actual internal
         * simulation domain nx-3 x ny-3 cells.
         */
        unsigned int nx, ny;

        /**
         * Time associated with this iteration/timestep.
         */
        float time;
};

#endif /* DATATYPES_H_ */
