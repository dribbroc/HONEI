/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * KPTypes.h
 *
 *  Created on: Mar 26, 2010
 *      Author: babrodtk
 */

#ifndef KPTYPES_H_
#define KPTYPES_H_

#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif
#include <iostream>
#include <cassert>





/**
 * This class describes the boundary conditions used in the simulation.
 */
class KPBoundaryConditions {
    public:
        enum TYPE {
            NONE=0,            //!< No boundary conditions
            WALL=1,            //!< Wall
            FIXED_DEPTH=2,     //!< Fixed water depth
            FIXED_DISCHARGE=3, //!< Fixed discharge (h*u / h*v)
            OPEN=4             //!< Open(? = nonreflective?) boundary conditions
        }; //!< This enum defines the different boundary conditions we can have.

    public:
        /**
         * Constructor that sets up boundary condition description
         * @param north_ North boundary condition
         * @param south_
         * @param east_
         * @param west_
         * @param north_arg_ Argument to BC. Depends on the type of condition what it is used for.
         * @param south_arg_
         * @param east_arg_
         * @param west_arg_
         * @return
         */
        KPBoundaryConditions(TYPE north_=NONE, TYPE south_=NONE,
                TYPE east_=NONE, TYPE west_=NONE,
                float north_arg_=0.0f, float south_arg_=0.0f,
                float east_arg_=0.0f, float west_arg_=0.0f) :
            north(north_), south(south_), east(east_), west(west_),
            north_arg(north_arg_), south_arg(south_arg_), east_arg(east_arg_), west_arg(west_arg_) {};

        /**
         * @return Boundary condition type for north boundary
         */
        inline const TYPE& getNorth() const { return north; }
        inline const TYPE& getSouth() const { return south; }
        inline const TYPE& getEast() const { return east; }
        inline const TYPE& getWest() const { return west; }

        /**
         * @return Argument to be used with boundary condition
         */
        inline const float& getNorthArg() const { return north_arg; }
        inline const float& getSouthArg() const { return south_arg; }
        inline const float& getEastArg() const { return east_arg; }
        inline const float& getWestArg() const { return west_arg; }

    private:
        TYPE north, south, east, west;                  //!< Boundary condition type for each interface
        float north_arg, south_arg, east_arg, west_arg; //!< Arguments for boundary conditions
};





/**
 * Helper function to easily print out boundary conditions
 */
inline std::ostream& operator<<(std::ostream& out, const KPBoundaryConditions::TYPE& type) {
    switch (type) {
        case KPBoundaryConditions::NONE:            out << "none"; break;
        case KPBoundaryConditions::WALL:            out << "wall"; break;
        case KPBoundaryConditions::FIXED_DEPTH:     out << "fixed depth"; break;
        case KPBoundaryConditions::FIXED_DISCHARGE: out << "fixed discharge"; break;
        case KPBoundaryConditions::OPEN:            out << "open"; break;
        default:                 out << "unknown"; break;
    }
    return out;
}




/**
 * Helper function to easily print out boundary conditions
 */
inline std::ostream& operator<<(std::ostream& out, const KPBoundaryConditions& bc) {
    out << "[n=" << bc.getNorth()
        << ", s=" << bc.getSouth()
        << ", e=" << bc.getEast()
        << ", w=" << bc.getWest() << "], ";
    out << "arg=[" << bc.getNorthArg() << ", "
        << bc.getSouthArg() << ", "
        << bc.getEastArg() << ", "
        << bc.getWestArg() << "]";
    return out;
}





/**
 * Struct that describes initial conditions we need to set up a simulation.
 * This object owns no data, and only keeps pointers to data, and other simulation
 * initial conditions.
 * FIXME: Document this struct thoroughly
 */
class KPInitialConditions {
    public:
        /**
         * Definition of different data locations. This is due to the use
         * of staggered grids in the Kurganov-Petrova scheme.
         */
        enum DATA_LOCATION {
            /**
             * Used to specify that data is to be interpreted as corner values.
             * This requires that the variable being describes has nx+1 by ny+1 values
             * as we have nx+1 by ny+1 intersections for nx by ny cells.
             */
            GRID_CELL_INTERSECTION,
            GRID_CELL_CENTER,        //!< Used to specify that data is to be interpreted as cell centered values
            GRID_CELL_UNKNOWN        //!< Unknown. This is an error typically
        };

        /**
         * Definition of different time ODE integrators
         */
        enum TIME_INTEGRATOR {
            EULER=0,            //!< Euler integration
            RUNGE_KUTTA_2=1,    //!< Second order accurate Runge-Kutta integration
            INTEGRATOR_UNKNOWN, //!< INTEGRATOR_UNKNOWN
        };

    public:

        /**
         * Constructor that enables us to have compile-time empty initial conditions.
         */
        KPInitialConditions() {
            reset();
        }

        virtual ~KPInitialConditions() { }

        /**
         * Constructor that sets all required variables for a simulation.
         * @param nx Number of grid cells in x direction
         * @param ny Number of grid cells in y direction
         * @param dx Grid cell spacing
         * @param dy Grid cell spacing
         * @param B_ Pointer to B-data
         * @param B_location_ Data location of B-data.
         * 		Note that only GRID_CELL_INTERSECTION is valid. This argument is included
         * 		so that users of the interface are aware that you have to take it into account.
         * @param U1_ Pointer to w-data
         * @param U2_ Pointer to hu-data
         * @param U3_ Pointer to hv-data
         * @param U_location Data location of U-data. Can be given at either centers or at intersections.
         * 		However, note that when given at intersections, it is assumed to be piecewise bi-linear,
         * 		and the value at grid cell centers is computed as the average of its four corners.
         */
        KPInitialConditions(unsigned int nx_, unsigned int ny_,
                float dx_, float dy_,
                KPBoundaryConditions bc_, TIME_INTEGRATOR ode_,
                DATA_LOCATION B_location_, float* B_,
                DATA_LOCATION U_location_, float* U1_, float* U2_=NULL, float* U3_=NULL) {
            reset();
            nx = nx_;
            ny = ny_;
            dx = dx_;
            dy = dy_;
            bc = bc_;
            ode = ode_;
            B_location = B_location_;
            B = B_;
            U_location = U_location_;
            U1 = U1_;
            U2 = U2_;
            U3 = U3_;
            h_desingularization_eps = 2.5e-3f*std::min(dx, dy);
        }

        /**
         * This function resets all variables to standard values.
         */
        inline void reset() {
            nx = 0;
            ny = 0;
            dx = 0.0f;
            dy = 0.0f;
            bc = KPBoundaryConditions();
            ode = INTEGRATOR_UNKNOWN;
            dt_scale = 1.0f;
            B_location = GRID_CELL_UNKNOWN;
            B = NULL;
            U_location = GRID_CELL_UNKNOWN;
            U1 = NULL;
            U2 = NULL;
            U3 = NULL;
            n = 0.0f;
            g = 9.80665f; //Standard gravity is not 9.81, but this value:)
            h_desingularization_eps = 0.0f;
        }

        /**
         * Checks whether these initialconditions are valid or not.
         */
        inline bool isValid() const {
            return (nx != 0
                    && ny != 0
                    && dx > 0.0f
                    && dy > 0.0f
                    && B != NULL
                    && U1 != NULL
                    && dt_scale > 0.0f
                    && B_location == GRID_CELL_INTERSECTION
                    && U_location != GRID_CELL_UNKNOWN
                    && ode != INTEGRATOR_UNKNOWN
                    && n >= 0.0f
                    && g > 0.0f
                    && h_desingularization_eps > 0.0f);
        }

        /**
         * @return Pointer to the bathymetry
         */
        float* getB() const {
            assert(isValid());
            return B;
        }

        /**
         * @return Pointer to water elevation
         */
        float* getU1() const {
            assert(isValid());
            return U1;
        }

        /**
         * @return Pointer to water discharge in u-direction
         */
        float* getU2() const {
            assert(isValid());
            return U2;
        }

        /**
         * @return Pointer to water discharge in v-direction
         */
        float* getU3() const {
            assert(isValid());
            return U3;
        }

        /**
         * @return Interpretation of bathymetry values. Should always return intersections
         */
        DATA_LOCATION getBLocation() const {
            assert(isValid());
            return B_location;
        }

        /**
         * @return Interpretation of [w, hu, hv] values
         */
        DATA_LOCATION getULocation() const {
            assert(isValid());
            return U_location;
        }

        /**
         * @return The time integrator type to use
         */
        TIME_INTEGRATOR getTimeIntegrator() const {
            assert(isValid());
            return ode;
        }

        /**
         * @return Number of values for the bathymetry pointer. This should always return
         * (nx+1) for a nx wide domain.
         */
        unsigned int getBNx() const {
            assert(isValid());
            return nx+1;
        }


        /**
         * @return Number of values for the bathymetry pointer. This should always return
         * (ny+1) for a ny wide domain.
         */
        unsigned int getBNy() const {
            assert(isValid());
            return ny+1;
        }

        /**
         * @return Number of values for the U-pointers
         */
        unsigned int getUNx() const {
            assert(isValid());
            switch (U_location) {
                case GRID_CELL_INTERSECTION: return nx+1;
                case GRID_CELL_CENTER: return nx;
                default: return 0;
            }
        }

        /**
         * @return Number of values for the U-pointers
         */
        unsigned int getUNy() const {
            assert(isValid());
            switch (U_location) {
                case GRID_CELL_INTERSECTION: return ny+1;
                case GRID_CELL_CENTER: return ny;
                default: return 0;
            }
        }

        /**
         * @return Width of domain in number of cells
         */
        const unsigned int& getNx() const {
            assert(isValid());
            return nx;
        }


        /**
         * @return Height of domain in number of cells
         */
        const unsigned int& getNy() const {
            assert(isValid());
            return ny;
        }

        /**
         * @return Grid cell spacing in x direction
         */
        const float& getDx() const {
            assert(isValid());
            return dx;
        }


        /**
         * @return Grid cell spacing in y direction
         */
        const float& getDy() const {
            assert(isValid());
            return dy;
        }

        /**
         * @return Scaling factor of Dt to ensure stability
         */
        const float& getDtScale() const {
            assert(isValid());
            return dt_scale;
        }

        /**
         * @return Desingularization parameter in the Kurganov-Petrova scheme.
         */
        const float& getDesingularizationEps() const {
            return h_desingularization_eps;
        }

        /**
         * @return Gravitational constant
         */
        const float& getG() const {
            assert(isValid());
            return g;
        }

        /**
         * @return Manning friction coefficient
         */
        const float& getN() const {
            assert(isValid());
            return n;
        }                  //!< Manning friction coefficient

        /**
         * @return Boundary conditions
         * @see KPBoundaryConditions
         */
        const KPBoundaryConditions& getBC() const {
            assert(isValid());
            return bc;
        } //!< Boundary conditions

        /**
         * Set grid cell spacing in x direction
         */
        void setDx(const float& dx_) {
            dx = dx_;
        }

        /**
         * Set grid cell spacing in y direction
         */
        void setDy(const float& dy_) {
            dy = dy_;
        }

        /**
         * Set gravitational constant
         */
        void setG(const float& g_) {
            g = g_;
        }

        /**
         * Set manning coefficient
         */
        void setN(const float& n_) {
            n = n_;
        }

        /**
         * Set boundary conditions
         * @see KPBoundaryConditions
         */
        void setBC(const KPBoundaryConditions& bc_) {
            bc = bc_;
        }

        /**
         * Set scaling of Dt. This is used to ensure that the simulation does not
         * become unstable. For 0 friction coefficient, this should typically be set
         * to a value less than 1.
         * @param scale
         */
        void setDtScale(const float& scale) {
            dt_scale = scale;
        }

        /**
         * Set desingularization epsilon used in the Kurganov-Petrova scheme. Because
         * fluxes use the value u = hu / h, we get extremely large u for very small h, and
         * the rounding error is enormous. This is fixed in the scheme by desingularizing the
         * hu / h calculation for values less than eps. eps is typically given in meters, so that
         * setting it to 1.0 will desingularize fluxes where the water depth is less than one meter.
         * @param eps
         */
        void setDesingularizationEps(const float& eps) {
            h_desingularization_eps = eps;
        }

    private:
        float* B;                 //!< Bathymetry
        float* U1;                //!< Water elevation
        float* U2;                //!< hu discharge
        float* U3;                //!< hv discharge
        DATA_LOCATION B_location; //!< How to interpret data
        DATA_LOCATION U_location; //!< How to interpret data
        TIME_INTEGRATOR ode;      //!< Time integrator

        KPBoundaryConditions bc;     //!< Boundary conditions
        unsigned int nx, ny;         //!< Number of initial grid cell intersections for bathymetry
        float dx, dy;                //!< Grid cell spacing
        float dt_scale;              //!< Scaling of dt to ensure stability
        float h_desingularization_eps; //!< Desingularization epsilon used in flux calculations

        float g; //!< Gravitational constant.
        float n; //!< Manning friction coefficient.
};





/**
 * Helper function to easily print out time integrator
 */
inline std::ostream& operator<<(std::ostream& out, const KPInitialConditions::TIME_INTEGRATOR& t) {
    switch(t) {
        case KPInitialConditions::EULER: out << "Euler"; break;
        case KPInitialConditions::RUNGE_KUTTA_2: out << "2nd order TVD Runge-Kutta"; break;
        default: out << "Unknown";
    }
    return out;
}




/**
 * Helper function to easily print out initial conditions
 */
inline std::ostream& operator<<(std::ostream& out, const KPInitialConditions& init) {
    out << "Domain size: [" << init.getNx() << "x" << init.getNy() << "], ";
    out << "cell size: [" << init.getDx() << "x" << init.getDy() << "]" << std::endl;
    out << "g = " << init.getG() << ", Mannings n = " << init.getN() << std::endl;
    out << "Boundary conditions: " << init.getBC() << std::endl;
    out << "Time integrator: " << init.getTimeIntegrator() << ", ";
    out << "Dt scaled by " << init.getDtScale();
    return out;
}



/**
 * Helper function to easily print out initialconditions
 */
inline std::ostream& operator<<(std::ostream& out, const KPInitialConditions* init) {
    out << *init;
    return out;
}






























#endif /* KPTYPES_H_ */
