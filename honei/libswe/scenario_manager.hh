/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBSWE_GUARD_SCENARIO_MANAGER_HH
#define LIBSWE_GUARD_SCENARIO_MANAGER_HH 1

#include <scenario.hh>
#include <honei/libutil/exception.hh>
#include <honei/libla/algorithm.hh>
/**
 * \file
 * Implementation of ScenarioManager and related (mostly Exception - based) classes.
 * \ingroup grplibswe
 **/
using namespace honei;
using namespace swe_solvers;

namespace honei
{

    /**
     * SWEError is the base class for all libswe errors.
     *
     *
     **/
    class SWEError : public Exception
    {
        protected:
            SWEError(const std::string & message) throw () : Exception(message)
        {
        }
    };

    /**
     * AllocationError is thrown, if user data location is unspecified.
     *
     **/
    class AllocationError : public SWEError
    {
        public:
            AllocationError(const std::string & p) throw () :SWEError(p + " has not been allocated! Use allocate_"+p+" !")
            {
            }
    };

    /**
     * EnvironmentalError is thrown, if environmental variables have not yet been set or have unexpected values.
     *
     **/
    class EnvironmentalError : public SWEError
    {
        public:
            EnvironmentalError(const std::string & p) throw () :SWEError(p)
            {
            }
    };

    /**
     * ScenarioManager has two basic functions: In the first place it controls allocation and setting of user data, secondly it provides functionality to load and save Scenarios using the HDF5 standard.
     *
     **/
    template<typename DataType_, typename SolverType_, typename BoundaryType_>
    class ScenarioManager
    {
    };

    /**
     * Explicit specialization for RelaxSolver.
     *
     **/
    template<typename DataType_>
    class ScenarioManager<DataType_, swe_solvers::RELAX, boundaries::REFLECT>
    {
        private:

            ///Needed scenario information:
            unsigned long _relax_entries;
            unsigned long _grid_width;
            unsigned long _grid_height;

            ///Flags for validation:
            bool _scenario_valid_flag;
            bool _allocated_scenario_flag;
            bool _allocated_scalarfields_flag;
            bool _allocated_relax_vectors_flag;
            bool _allocated_bottom_slopes_flag;
            bool _environment_delta_set_flag;
            bool _environment_manning_set_flag;
            bool _environment_epsilon_set_flag;

            ///Our scenario
            Scenario<DataType_, swe_solvers::RELAX, boundaries::REFLECT> * _current_scenario;

            /**
             * Sets all needed discretization environmental variables.
             *
             * \param dx The stepsize in x direction.
             * \param dy The stepsize in y direction.
             * \param dt The size of the timestep.
             **/
            void _set_environment_delta(DataType_ dx, DataType_ dy, DataType_ dt)
            {
                this->_current_scenario->delta_x = dx;
                this->_current_scenario->delta_y = dy;
                this->_current_scenario->delta_t = dt;

                this->_environment_delta_set_flag = true;
            }

            /**
             * Sets the Manning constant.
             *
             * \param mn The Manning constant to be used.
             *
             **/
            void _set_environment_manning(DataType_ mn)
            {
                this->_current_scenario->manning_n = mn;
                this->_environment_manning_set_flag = true;
            }

            /**
             * Sets epsilon.
             *
             * \param eps Our epsilon relaxation parameter.
             **/
            void _set_environment_epsilon(double eps)
            {
                this->_current_scenario->eps = eps;
                this->_environment_epsilon_set_flag = true;
            }

        public:

            ScenarioManager() :
                _scenario_valid_flag(false),
                _allocated_scenario_flag(false),
                _allocated_scalarfields_flag(false),
                _allocated_relax_vectors_flag(false),
                _allocated_bottom_slopes_flag(false),

                _environment_delta_set_flag(false),
                _environment_manning_set_flag(false),
                _environment_epsilon_set_flag(false)
            {
            }

            /**
             * Assigns the scenario.
             *
             * \param scenario Pointer to a Scenario object.
             **/
            void allocate_scenario(Scenario<DataType_, swe_solvers::RELAX, boundaries::REFLECT> * scenario)
            {
                this->_current_scenario = scenario;
                this->_allocated_scenario_flag = true;
            }

            /**
             * Converts all data from one precision to another - this will work only with target scenarios members all allocated, due to performance reasons, this wont be checked and allocation cant be done here!
             *
             * \param target Reference to target scenario object.
             * \param source Reference to source scenario object.
             **/
            template <typename TargetDT_, typename SourceDT_>
            static void convert_scenario(Scenario<TargetDT_, swe_solvers::RELAX, boundaries::REFLECT> & target,
                                         Scenario<SourceDT_, swe_solvers::RELAX, boundaries::REFLECT> & source)
            {
                convert(*(target.height), *(source.height));
                convert(*(target.bottom), *(source.bottom));
                convert(*(target.x_veloc), *(source.x_veloc));
                convert(*(target.y_veloc), *(source.y_veloc));
                convert(*(target.bottom_slopes_x), *(source.bottom_slopes_x));
                convert(*(target.bottom_slopes_y), *(source.bottom_slopes_y));

                convert(*(target.u), *(source.u));
                convert(*(target.v), *(source.v));
                convert(*(target.w), *(source.w));
                convert(*(target.c), *(source.c));
                convert(*(target.d), *(source.c));

                target.delta_x = (TargetDT_)source.delta_x;
                target.delta_y = (TargetDT_)source.delta_y;
                target.delta_t = (TargetDT_)source.delta_t;
                target.manning_n = (TargetDT_)source.manning_n;
            }

            /**
             * Assigns all scalarfields.
             *
             * \param h Pointer to our height field.
             * \param b Pointer to our bottom field.
             * \param u Pointer to our velocity field (x-direction)
             * \param v Pointer to our velocity field (y-direction)
             *
             **/
            void allocate_scalarfields(DenseMatrix<DataType_> * h, DenseMatrix<DataType_> * b, DenseMatrix<DataType_> * u, DenseMatrix<DataType_> * v)
            {
                this->_current_scenario->height = h;
                this->_grid_width = h->columns();
                this->_grid_height = h->rows();
                this->_relax_entries = 3*((_grid_width * _grid_height)+4*(_grid_width +_grid_height+4));

                if(h->rows() != b->rows())
                    throw MatrixRowsDoNotMatch(b->rows(), _grid_height);
                if(h->columns() != b->columns())
                    throw MatrixColumnsDoNotMatch(b->columns(), _grid_width);
                this->_current_scenario->bottom = b;

                if(h->rows() != u->rows())
                    throw MatrixRowsDoNotMatch(u->rows(), _grid_height);
                if(h->columns() != u->columns())
                    throw MatrixColumnsDoNotMatch(u->columns(), _grid_width);
                this->_current_scenario->x_veloc = u;

                if(h->rows() != v->rows())
                    throw MatrixRowsDoNotMatch(v->rows(), _grid_height);
                if(h->columns() != v->columns())
                    throw MatrixColumnsDoNotMatch(v->columns(), _grid_width);

                this->_current_scenario->y_veloc = v;
                this->_allocated_scalarfields_flag = true;
            }

            /**
             * Assigns all relax vectors needed by the solver.
             *
             * \param u Pointer to relaxation vector u.
             * \param v Pointer to relaxation vector v.
             * \param w Pointer to relaxation vector w.
             * \param c Pointer to relaxation vector c.
             * \param d Pointer to relaxation vector d.
             *
             **/
            void allocate_relax_vectors(DenseVector<DataType_> * u, DenseVector<DataType_> * v, DenseVector<DataType_> * w, DenseVector<DataType_> * c, DenseVector<DataType_> * d)
            {
                if(u->size() != _relax_entries)
                    throw VectorSizeDoesNotMatch(u->size(), _relax_entries);
                this->_current_scenario->u = u;
                if(v->size() != _relax_entries)
                    throw VectorSizeDoesNotMatch(v->size(), _relax_entries);
                this->_current_scenario->v = v;
                if(w->size() != _relax_entries)
                    throw VectorSizeDoesNotMatch(w->size(), _relax_entries);
                this->_current_scenario->w = w;
                if(c->size() != 3)
                    throw VectorSizeDoesNotMatch(c->size(), 3);
                this->_current_scenario->c = c;
                if(d->size() != 3)
                    throw VectorSizeDoesNotMatch(d->size(), 3);
                this->_current_scenario->d = d;
                this->_allocated_relax_vectors_flag = true;
            }

            /**
             * Assigns all slope vectors.
             *
             * \param bx Pointer to slope vector x.
             * \param by Pointer to slope vector y.
             *
             **/
            void allocate_bottom_slopes(DenseVector<DataType_> * bx, DenseVector<DataType_> * by)
            {
                if(bx->size() != _relax_entries/3)
                    throw VectorSizeDoesNotMatch(bx->size(), _relax_entries/3);
                this->_current_scenario->bottom_slopes_x = bx;
                if(by->size() != _relax_entries/3)
                    throw VectorSizeDoesNotMatch(by->size(), _relax_entries/3);
                this->_current_scenario->bottom_slopes_y = by;
                this->_allocated_bottom_slopes_flag = true;
            }

            /**
             * Sets all environmental variables - see private methods for details.
             **/
            void set_environmental_variables(DataType_ dx, DataType_ dy, DataType_ dt, DataType_ mn, DataType_ eps)
            {
                _set_environment_delta(dx, dy, dt);
                _set_environment_manning(mn);
                _set_environment_epsilon(eps);
            }

            /**
             * Validates the current scenario.
             *
             **/
            bool validate()
            {
                if(!_allocated_scenario_flag)
                    throw AllocationError("scenario");

                if(!_allocated_scalarfields_flag)
                    throw AllocationError("scalarfields");

                if(!_allocated_relax_vectors_flag)
                    throw AllocationError("relax_vectors");

                if(!_allocated_bottom_slopes_flag)
                    throw AllocationError("bottom_slopes");

                if(!_environment_delta_set_flag)
                    throw EnvironmentalError("Environment unset! Call to set_environmental_variables needed!");

                this-> _scenario_valid_flag = true;
                return true;
            }

            /**
             * Saves current scenario in HDF5 format.
             *
             * \param filename The name of the HDF5 file.
             *
             **/
            void save(const std::string filename)
            {
            }

            /**
             * Loads scenario from an existing HDF5 file.
             *
             * \param filename The name of the HDF5 file.
             *
             **/
            void load(const std::string filename)
            {
                //TODO
            }
    };
}
#endif
