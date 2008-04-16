/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
 *
 * This file is part of the Graph C library. LibGraph is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBGRAPH_GUARD_POSITIONS_HH
#define LIBGRAPH_GUARD_POSITIONS_HH 1

#include <honei/libgraph/node_distance.hh>
#include <honei/libgraph/position-impl.hh>
#include <honei/libgraph/breadth_first_search.hh>
#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/difference.hh>
#include <honei/libla/dot_product.hh>
#include <honei/libla/element_inverse.hh>
#include <honei/libla/element_product.hh>
#include <honei/libla/matrix_error.hh>
#include <honei/libla/product.hh>
#include <honei/libla/reduction.hh>
#include <honei/libla/sum.hh>
#include <honei/libla/scale.hh>
#include <honei/libla/vector.hh>
#include <honei/libgraph/abstract_graph.hh>
#include <honei/libgraph/graph_error.hh>
#include <honei/libutil/tags.hh>

#include <cmath>
#include <tr1/memory>
#include <iostream>
#include <fstream> 
#include <queue>
/**
 * \file
 *
 * Interface declarations for position calculation classes.
 *
 * \ingroup grplibgraph
 **/
 namespace honei
{
    namespace methods
    {
        /**
         * Implementation class template for the varying methods.
         */
        template <typename Tag_, typename DataType_, typename GraphTag_>
        class Implementation;

    }

        /**
         * class template for statistics.
         */
    template <typename Tag_> class Statistics
    {
        public:
            template <typename DataType_>
            static void value(const unsigned long _number_of_iterations,
            const DenseMatrix<DataType_> & square_dist, const DenseMatrix<DataType_> & _weights_of_edges, const DataType_ & max_force, 
            DenseVector<DataType_> & statistic_result)
            {
                CONTEXT("When calculating the average edge length and the standard deviation of lengths:");

                // statistic_result containing number of iterations, average length and standard deviation
                statistic_result[0] = _number_of_iterations;
                statistic_result[1] = 0;
                statistic_result[2] = 0;
                statistic_result[3] = max_force;
                unsigned long number_of_edges(0);

                // calculating average edge length
                for (unsigned long i(0); i < square_dist.rows(); ++i)
                {
                    for (unsigned long j(i+1); j < square_dist.columns(); ++j)
                    {
                        if(_weights_of_edges[i][j] > 0)
                        {
                            statistic_result[1] += sqrt(square_dist[i][j]);
                            number_of_edges++;
                        }
                    }
                }
                statistic_result[1] /= number_of_edges;

                // calculating standard deviation of lengths
                for (unsigned long i(0); i < square_dist.rows(); ++i)
                {
                    for (unsigned long j(i+1); j < square_dist.columns(); ++j)
                    {
                        if(_weights_of_edges[i][j] > 0)
                        {
                            DataType_ diff(sqrt(square_dist[i][j]) - statistic_result[1]);
                            statistic_result[2] += diff*diff;
                        }
                    }
                }
                statistic_result[2] = sqrt(statistic_result[2] / (number_of_edges - 1));
            }

            template <typename DataType_>
            static void value(const unsigned long _number_of_iterations,
            const SparseMatrix<DataType_> & _weights_of_edges, const DataType_ & max_force, const DenseMatrix<DataType_> & _coordinates, 
            DenseVector<DataType_> & statistic_result)
            {
                CONTEXT("When calculating the average edge length and the standard deviation of lengths:");

                // statistic_result containing number of iterations, average length and standard deviation
                statistic_result[0] = _number_of_iterations;
                statistic_result[1] = 0;
                statistic_result[2] = 0;
                statistic_result[3] = max_force;
                unsigned long number_of_edges(0);

                // calculating average edge length
                for (unsigned long i(0); i < _weights_of_edges.rows(); ++i)
                {
                    for (unsigned long j(i+1); j < _weights_of_edges.columns(); ++j)
                    {
                        if(_weights_of_edges[i][j] > 0)
                        {
                            DenseVector<DataType_> v(_coordinates.columns());
                            TypeTraits<DataType_>::copy(_coordinates[i].elements(), v.elements(), _coordinates.columns());
                            DataType_ d(Norm<vnt_l_two, true, tags::CPU>::value(Difference<tags::CPU>::value(v, _coordinates[j])));
                            statistic_result[1] += d;
                            number_of_edges++;
                        }
                    }
                }
                statistic_result[1] /= number_of_edges;

                // calculating standard deviation of lengths
                for (unsigned long i(0); i < _weights_of_edges.rows(); ++i)
                {
                    for (unsigned long j(i+1); j < _weights_of_edges.columns(); ++j)
                    {
                        if(_weights_of_edges[i][j] > 0)
                        {
                            DenseVector<DataType_> v(_coordinates.columns());
                            TypeTraits<DataType_>::copy(_coordinates[i].elements(), v.elements(), _coordinates.columns());
                            DataType_ d(Norm<vnt_l_two, true, tags::CPU>::value(Difference<tags::CPU>::value(v, _coordinates[j])));
                            DataType_ diff(d - statistic_result[1]);
                            statistic_result[2] += diff*diff;
                        }
                    }
                }
                statistic_result[2] = sqrt(statistic_result[2] / (number_of_edges - 1));
            }
    };

    template <typename Tag_ , typename DataType_, typename GraphTag_> class Positions
    {
        private:
            /// Our implementation pattern.
            methods::Implementation<Tag_ ,DataType_, GraphTag_> * _imp;

            inline void elements_check(const SparseMatrix<DataType_> & weights_of_edges, bool & trace, bool & symmetry, bool & positiv)
            {
                for (typename Matrix<DataType_>::ConstElementIterator e(weights_of_edges.begin_non_zero_elements()),
                        e_end(weights_of_edges.end_non_zero_elements()); e != e_end ; ++e)
                {
                    if (*e != weights_of_edges(e.column(), e.row())) symmetry = false;
                    if (e.column() == e.row()) trace = false;
                    if (*e < 0) positiv = false;
                }
            }

            inline void elements_check(const DenseVector<DataType_> & weights_of_nodes, bool & positiv)
            {
                for (typename Vector<DataType_>::ConstElementIterator e(weights_of_nodes.begin_elements()),
                e_end(weights_of_nodes.end_elements()); e != e_end ; ++e)
                {
                    if (*e < 0) positiv = false;
                }
            }

        public:  
            Positions(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                _imp(new methods::Implementation<Tag_, DataType_, GraphTag_>(coordinates, weights_of_nodes, weights_of_edges))
            {
                bool symmetric(true);
                bool trace(true);
                bool positiv(true);
                bool positiv_2(true);
                if (coordinates.rows() != weights_of_edges.columns())
                    throw MatrixColumnsDoNotMatch(weights_of_edges.columns(), coordinates.columns());

                if (! weights_of_edges.square())
                    throw MatrixIsNotSquare(weights_of_edges.rows(), weights_of_edges.columns());

                if (weights_of_nodes.size() != coordinates.rows())
                    throw VectorSizeDoesNotMatch(weights_of_nodes.size(), coordinates.columns());

                elements_check(weights_of_edges, trace, symmetric, positiv);
                if (! symmetric)
                    throw GraphError("Matrix of edge-weights must be symmetric");

                if (! trace)
                    throw GraphError("Matrix of edge-weights has non vanishing trace");

                if (! positiv)
                    throw GraphError("Elements in edge-weights matrix must be positiv");

                elements_check(weights_of_nodes, positiv_2);
                if (! positiv_2)
                    throw GraphError("Elements in node-weights vector must be positiv");
            }

            Positions(AbstractGraph<DataType_> & graph, DataType_ edge_length) :
                _imp(new methods::Implementation<Tag_, DataType_, GraphTag_>(graph, edge_length))
            {
                /* Following lines not needed: Graph assembles matrices correctly - I hope so, at least.
                 *
                if (graph.coordinates->columns() != graph.cweights_of_edges.columns())
                    throw MatrixColumnsDoNotMatch(weights_of_edges.columns(), coordinates.columns());

                if (! weights_of_edges.square())
                    throw MatrixIsNotSquare(weights_of_edges.rows(), weights_of_edges.columns());

                if (weights_of_nodes.size() != coordinates.columns())
                    throw VectorSizeDoesNotMatch(weights_of_nodes.size(), coordinates.columns());
                    */
#if 0
                \todo Implement Trace
                if (Trace<>::value(weights_of_edges) > std::numeric_limits<DataType_>::epsilon())
                    throw MatrixHasNonVanishingTrace;
#endif
            }

            ~Positions()
            {
                delete _imp;
            }

            void update(DataType_ eps, unsigned long timeout)
            {
                _imp->init();
                do
                {
                    --timeout;
                }
                while ((timeout > 0) && (_imp->value(eps) > eps));
            }

            void init()
            {
                _imp->init();
            }

            DataType_ step()
            {
                DataType_ dt(0);
                return _imp->value(dt);
            }

            DenseMatrix<DataType_> & coordinates()
            {
                return _imp->_coordinates;
            }

            DataType_ max_node_force()
            {
                return _imp->_max_force;
            }

            unsigned long number_of_iterations()
            {
                return _imp->_number_of_iterations;
            }
            
            inline  DenseVector<DataType_> step_width()
            {
                return _imp->step_width();
            }

            inline void step_width(DataType_ size)
            {
                _imp->step_width(size);
            }

            inline void step_width_factors(DataType_ factor_1, DataType_ factor_2)
            {
                _imp->step_width_factors(factor_1, factor_2);
            }

            void noise_duration(DataType_ size)
            {
                _imp->noise_duration(size);
            }

            void statistic(unsigned long size)
            {
                _imp->statistic(size);
            }

            std::queue< DenseVector<DataType_> > statistic()
            {
                return _imp->statistic();
            }

            void get_statistic()
            {
                _imp->get_statistic();
            }
            
            void get_statistic(char filename[])
            {
                _imp->get_statistic(filename);
            }
    };

    namespace methods
    {
        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, WeightedFruchtermanReingold>
        {
            private:
                ///  position matrix - the coordinates of each node
                DenseMatrix<DataType_> _coordinates;

                ///  weight vector - the weights of each node
                const DenseVector<DataType_> _weights_of_nodes;

                ///  edge weight matrix - the weights of each edge
                const DenseMatrix<DataType_> _weights_of_edges;

                /// parameter matrix for repulsive forces
                DenseMatrix<DataType_> _repulsive_force_parameter;

                /// parameter matrix for attractive forces
                DenseMatrix<DataType_> _attractive_force_parameter;

                /// number of iterations
                unsigned long _number_of_iterations;

                /// the maximal force
                DataType_ _max_force;

                /// force direction
                DenseMatrix<DataType_> _force_direction;

                /// step width
                DenseVector<DataType_> _step_width;

                /// step width factors: factor_1 increases step_width (if node is moving in correct direction), factor_2 reduces step_width (if node is oscillating or rotating)
                DataType_ _step_width_factor_1, _step_width_factor_2;

                /// _noise_duration - how long is a noise added to the forces (related to _step_with)?
                DataType_ _noise_duration;

                /// the _repulsive_force_range defines the range of the repulsive force
                DataType_ _repulsive_force_range;

                /// statistic_queue contains statistic values (number of iterations, average edge length, standard deviation of lengths and the maximum node force)
                std::queue< DenseVector<DataType_> > statistic_queue;

                /// _statistic_step defines the iterations, in which statistic values will be calculated
                unsigned long _statistic_step;

            public:
                friend class Positions<Tag_, DataType_, WeightedFruchtermanReingold>;
                
                inline DenseVector<DataType_> step_width()
                {
                    return _step_width;
                }

                inline void step_width(DenseVector<DataType_> size)
                {
                    _step_width = size;
                }

                inline void step_width_factors(DataType_ factor_1, DataType_ factor_2)
                {
                    _step_width_factor_1 = factor_1;
                    _step_width_factor_2 = factor_2;
                }

                inline void noise_duration(DataType_ size)
                {
                    _noise_duration = size;
                }

                inline void statistic(unsigned long size)
                {
                    _statistic_step = size;
                }

                inline std::queue< DenseVector<DataType_> > statistic()
                {
                    return statistic_queue;
                }

                inline void get_statistic()
                {
                    std::cout<<"number of iterations: "<<"average edge length: "<<"standard deviation of lengths: "<<"maximum force: "<< std::endl;
                    while (!statistic_queue.empty())
                    {                        
                        DenseVector<DataType_> current_value(statistic_queue.front());
                        std::cout<<current_value[0]<<" "<<current_value[1]<<" "<<current_value[2]<<" "<<current_value[3]<< std::endl;
                        statistic_queue.pop();
                    }
                }
                
                inline void get_statistic(char filename[])
                {
                    std::ofstream file(filename); 
                    file<<"number_of_iterations: "<<"average_edge_length: "<<"standard_deviation_of_lengths: "<<"maximum_force: "<< std::endl;
                    while (!statistic_queue.empty())
                    {                        
                        DenseVector<DataType_> current_value(statistic_queue.front());
                        file<<current_value[0]<<" "<<current_value[1]<<" "<<current_value[2]<<" "<<current_value[3]<< std::endl;
                        statistic_queue.pop();
                    }
                }

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _repulsive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), DataType_(0)),
                    _attractive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), DataType_(0)),
                    _number_of_iterations(1),
                    _max_force(0),
                    _force_direction(coordinates.rows(), coordinates.columns(), DataType_(0)),
                    _step_width(weights_of_nodes.size()),
                    _repulsive_force_range(DataType_(0)),
                    _noise_duration(DataType_(0)),
                    _statistic_step(0),
                    _step_width_factor_1(1.5),
                    _step_width_factor_2(0.5)
                {
                    // Calculate parameter matrix for repulsive forces and attractive forces, _repulsive_force_range and _step_width
                    DataType_ _max_ideal_length(0);
                    typename MutableMatrix<DataType_>::ElementIterator f(_attractive_force_parameter.begin_elements());
                    typename Matrix<DataType_>::ConstElementIterator g(_weights_of_edges.begin_elements());
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_repulsive_force_parameter.begin_elements()),
                        e_end(_repulsive_force_parameter.end_elements()); e != e_end ; ++e, ++f, ++g)
                    {
                        DataType_ node_weight(_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]);
                        if (*g > std::numeric_limits<DataType_>::epsilon())
                        {
                            DataType_ length_of_edge(sqrt(node_weight) / *g);
                            if ((length_of_edge) > _max_ideal_length) _max_ideal_length = length_of_edge;
                        }
                        *e = node_weight * node_weight;
                        *f = *g * *g * *g * *g;
                    }
                    _noise_duration = _max_ideal_length / 40;
                    _repulsive_force_range = _max_ideal_length * _weights_of_edges.rows();
                    DataType_ s_w(_repulsive_force_range / 20);

                    for (typename Vector<DataType_>::ElementIterator e(_step_width.begin_elements()),
                        e_end(_step_width.end_elements()); e != e_end ; ++e)
                        {
                            *e = s_w;
                        }
                }

                Implementation(AbstractGraph<DataType_> & graph, DataType_ edgeLength) :
                    _coordinates(*graph.coordinates()),
                    _weights_of_nodes(*graph.node_weights()),
                    _weights_of_edges(*graph.edges()),
                    _repulsive_force_parameter(_weights_of_edges.columns(), _weights_of_edges.rows(), DataType_(0)),
                    _attractive_force_parameter(_weights_of_edges.columns(), _weights_of_edges.rows(), DataType_(0)),
                    _number_of_iterations(1),
                    _max_force(0),
                    _force_direction(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _step_width(_weights_of_nodes.size()),
                    _repulsive_force_range(DataType_(0)),
                    _noise_duration(DataType_(0)),
                    _statistic_step(0),
                    _step_width_factor_1(1.5),
                    _step_width_factor_2(0.5)
                {
                    // Calculate parameter matrix for repulsive forces and attractive forces, _repulsive_force_range and _step_width
                    DataType_ _max_ideal_length(0);
                    typename MutableMatrix<DataType_>::ElementIterator f(_attractive_force_parameter.begin_elements());
                    typename Matrix<DataType_>::ConstElementIterator g(_weights_of_edges.begin_elements());
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_repulsive_force_parameter.begin_elements()),
                        e_end(_repulsive_force_parameter.end_elements()); e != e_end ; ++e, ++f, ++g)
                    {
                        DataType_ node_weight(_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]);
                        if (*g > std::numeric_limits<DataType_>::epsilon())
                        {
                            DataType_ length_of_edge(sqrt(node_weight) / *g);
                            if ((length_of_edge) > _max_ideal_length) _max_ideal_length = length_of_edge;
                        }
                        *e = (graph.same_timeslice(e.row(), e.column())) ?
                           node_weight * node_weight: 0;
                        *f = *g * *g * *g * *g;
                    }
                    _noise_duration = _max_ideal_length / 40;
                    _repulsive_force_range = _max_ideal_length * _weights_of_edges.rows();
                    DataType_ s_w(_repulsive_force_range / 20);

                    for (typename Vector<DataType_>::ElementIterator e(_step_width.begin_elements()),
                        e_end(_step_width.end_elements()); e != e_end ; ++e)
                        {
                            *e = s_w;
                        }
                }

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2) and _square_dist = (d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0));
                    DenseMatrix<DataType_> square_dist(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0));
                    NodeDistance<Tag_>::value(_coordinates, _weights_of_edges, square_dist, inv_square_dist, _repulsive_force_range);
                    DenseMatrix<DataType_> square_dist_copy(square_dist.copy());

                    // Calculate inv_square_dist = Mul(inv_square_dist,  _repulsive_force_parameter)
                    ElementProduct<Tag_>::value(inv_square_dist, _repulsive_force_parameter);

                    // Calculate square_dist = Mul(square_dist,  _attractive_force_parameter)
                    ElementProduct<Tag_>::value(square_dist, _attractive_force_parameter);

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> sum_vec(Reduction<rt_sum, Tag_>::value(square_dist));

                    // Calculate the same stuff for inv_square_dist
                    DenseVector<DataType_> sum_vec_inv(Reduction<rt_sum, Tag_>::value(inv_square_dist));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec.size(), sum_vec);
                    Difference<Tag_>::value(diag, square_dist);
                    Scale<Tag_>::value(square_dist, DataType_(-1));

                    // Calculating attractive_forces = (square_dist * _coordinates)
                    DenseMatrix<DataType_> attractive_forces(Product<Tag_>::value(square_dist, _coordinates));

                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<Tag_>::value(inv_square_dist, DataType_(-1));

                    // Calculating diff = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv.size(), sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<Tag_>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (diff * _coordinates)
                    DenseMatrix<DataType_> repulsive_forces(Product<Tag_>::value(diff, _coordinates));

                    // Calculate the maximal force and the result forces
                    DataType_ result(0);
                    Sum<Tag_>::value(attractive_forces, repulsive_forces);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(attractive_forces[i.index()]);
                        result = std::max(result, *i);
                    }

                    // Calculate the new _step_width
                    DenseMatrix<DataType_> scaled_forces(_force_direction.rows(), _force_direction.columns(), DataType_(0));
                    for (typename MutableMatrix<DataType_>::ElementIterator e(scaled_forces.begin_elements()),
                        e_end(scaled_forces.end_elements()), k(attractive_forces.begin_elements()); e != e_end ; ++e, ++k)
                    {
                        resulting_forces[e.row()] > 0  ? *e = *k / resulting_forces[e.row()] :
                        *e = 0;
                    }
                    if (_number_of_iterations > 1)
                    {
                        for (typename Vector<DataType_>::ElementIterator e(_step_width.begin_elements()),
                                e_end(_step_width.end_elements()); e != e_end ; ++e)
                                {
                                        DataType_ prod( DotProduct<Tag_>::value(_force_direction[e.index()], scaled_forces[e.index()]) );
                                        if (prod > 0.8) *e *=_step_width_factor_1;
                                        if ( (prod < -0.8) || (fabs(prod) < 0.2) ) *e *=_step_width_factor_2;
                                }
                    }

                    // Calculate the new positions by using scaled forces
                    DataType_ noise(1);
                    DataType_ delta(0);
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_coordinates.begin_elements()),
                        e_end(_coordinates.end_elements()), k(scaled_forces.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        if ( e.column() == 0 )
                        {
                            _step_width[e.row()] > _noise_duration ? noise = (7.5 + ( rand() % 5 ) ) / 10 : noise = 1;
                            delta = std::min(_step_width[e.row()], resulting_forces[e.row()]);
                        }
                        result > eps ? *e = *e + delta * noise * *k : 0;
                    }

                    // Calculate statistic values
                    if ( (_statistic_step > 0) && (_number_of_iterations-1 == 0 || ((_number_of_iterations-1) % _statistic_step) == 0 ) )
                    {
                        DenseVector<DataType_> statistic_values(4, DataType_(0));
                        Statistics<Tag_>::value(_number_of_iterations-1, square_dist_copy, _weights_of_edges, result, statistic_values);
                        statistic_queue.push(statistic_values);
                    }

                    _force_direction = scaled_forces;
                    _number_of_iterations++;
                    _max_force = result;

                    return result;
                }

                void init()
                {
                        _number_of_iterations = 1;
                }
        };

        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, WeightedKamadaKawai>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> _coordinates;

                /// weight vector - the weights of each node
                const DenseVector<DataType_> _weights_of_nodes;

                /// edge weight matrix - the weights of each edge
                const SparseMatrix<DataType_> _weights_of_edges;

                /// graph distance matrix
                DenseMatrix<DataType_> _graph_distance;

                /// spring forces
                DenseMatrix<DataType_> _spring_forces;

                /// All forces of a node
                DenseMatrix<DataType_> node_forces;

                /// the maximal force
                DataType_ _max_force;

                /// number of iterations
                int _number_of_iterations;

                /// index of node with maximal force
                int _max_node;

                /// force direction
                DenseMatrix<DataType_> _force_direction;

                /// vector of step width
                DenseVector<DataType_> _step_width;

                /// step width factors: factor_1 increases step_width (if node is moving in correct direction), factor_2 reduces step_width (if node is oscillating or rotating)
                DataType_ _step_width_factor_1, _step_width_factor_2;

                /// _noise_duration - how long is a noise added to the forces (related to _step_with)?
                DataType_ _noise_duration;

                /// statistic_queue contains statistic values (number of iterations, average edge length, standard deviation of lengths and the maximum node force)
                std::queue< DenseVector<DataType_> > statistic_queue;

                /// _statistic_step defines the iterations, in which statistic values will be calculated
                unsigned long _statistic_step;

            public:
                friend class Positions<Tag_, DataType_, WeightedKamadaKawai>;

                inline DenseVector<DataType_> step_width()
                {
                    return _step_width;
                }

                inline void step_width(DenseVector<DataType_> size)
                {
                    _step_width = size;
                }

                inline void step_width_factors(DataType_ factor_1, DataType_ factor_2)
                {
                    _step_width_factor_1 = factor_1;
                    _step_width_factor_2 = factor_2;
                }

                inline void statistic(unsigned long size)
                {
                    _statistic_step = size;
                }

                inline std::queue< DenseVector<DataType_> > statistic()
                {
                    return statistic_queue;
                }

                inline void get_statistic()
                {
                    std::cout<<"number of iterations: "<<"average edge length: "<<"standard deviation of lengths: "<<"maximum force: "<< std::endl;
                    while (!statistic_queue.empty())
                    {                        
                        DenseVector<DataType_> current_value(statistic_queue.front());
                        std::cout<<current_value[0]<<" "<<current_value[1]<<" "<<current_value[2]<<" "<<current_value[3]<< std::endl;
                        statistic_queue.pop();
                    }
                }
                
                inline void get_statistic(char filename[])
                {
                    std::ofstream file(filename); 
                    file<<"number_of_iterations: "<<"average_edge_length: "<<"standard_deviation_of_lengths: "<<"maximum_force: "<< std::endl;
                    while (!statistic_queue.empty())
                    {                        
                        DenseVector<DataType_> current_value(statistic_queue.front());
                        file<<current_value[0]<<" "<<current_value[1]<<" "<<current_value[2]<<" "<<current_value[3]<< std::endl;
                        statistic_queue.pop();
                    }
                }

                inline void noise_duration(DataType_ size)
                {
                    _noise_duration = size;
                }

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _graph_distance(weights_of_edges.rows(), weights_of_edges.columns(), DataType_(0)),
                    _spring_forces(coordinates.rows(), coordinates.columns(), DataType_(0)),
                    _max_force(0),
                    _number_of_iterations(1),
                    _max_node(0),
                    node_forces(coordinates.rows() * coordinates.rows(), coordinates.columns(), DataType_(0)),
                    _step_width(coordinates.rows(), DataType_(0)),
                    _force_direction(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _statistic_step(0),
                    _noise_duration(0),
                    _step_width_factor_1(1.05),
                    _step_width_factor_2(0.35)
                {
                    // Using BFS to calculate the graph distance matrix
                    if (! BreadthFirstSearch<Tag_>::value(_graph_distance, _weights_of_nodes, _weights_of_edges))
                        throw GraphError("Graph must be coherent");
                }

                Implementation(AbstractGraph<DataType_> & graph, DataType_ edge_length) :
                    _coordinates(*graph.coordinates()),
                    _weights_of_nodes(*graph.node_weights()),
                    _weights_of_edges(*graph.edges()),
                    _graph_distance(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0)),
                    _spring_forces(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _max_force(0),
                    _number_of_iterations(1),
                    _max_node(0),
                    node_forces(_coordinates.rows() * _coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _step_width(_coordinates.rows(), DataType_(0)),
                    _force_direction(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _statistic_step(0),
                    _noise_duration(0),
                    _step_width_factor_1(1.05),
                    _step_width_factor_2(0.35)
                {
                    // Using BFS to calculate the graph distance matrix
                    if (! BreadthFirstSearch<Tag_>::value(_graph_distance, _weights_of_nodes, _weights_of_edges, graph))
                        throw GraphError("Graph must be coherent");
                }

                void init()
                {
                    _number_of_iterations = 1;
                    
                    // Calculate stepwidth
                    DataType_ stepwidth(0);
                    for (typename Matrix<DataType_>::ConstElementIterator e(_weights_of_edges.begin_non_zero_elements()),
                        e_end(_weights_of_edges.end_non_zero_elements()); e != e_end ; ++e)
                    {
                        DataType_ auxiliary(sqrt(_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]) / *e * _coordinates.rows() / 20);
                        if (auxiliary > stepwidth)
                        {
                            stepwidth = auxiliary;
                        }
                    }
                    _noise_duration = stepwidth / 20;

                    // Fill vector of _step_width with stepwidth
                    for (typename DenseVector<DataType_>::ElementIterator i(_step_width.begin_elements()), i_end(_step_width.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = stepwidth;
                    }

                    // Calculate square _graph_distance
                    ElementProduct<Tag_>::value(_graph_distance, _graph_distance);

                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<Tag_>::value(_coordinates));
                    DenseMatrix<DataType_> square_dist_copy(square_dist.copy());

                    // _spring_force_parameters = 1 / sum( _graph_distance, square_dist)
                    DenseMatrix<DataType_> _spring_force_parameters(_graph_distance.copy());
                    _spring_force_parameters = Sum<Tag_>::value(_spring_force_parameters, square_dist);
                    ElementInverse<Tag_>::value(_spring_force_parameters);

                    // Calculate square_dist = square_dist * 2
                    Scale<Tag_>::value(square_dist, DataType_(2));

                    // _spring_force_parameters = _spring_force_parameters * square_dist -1
                    ElementProduct<Tag_>::value(_spring_force_parameters, square_dist);
                    Sum<Tag_>::value(_spring_force_parameters, DataType_(-1));

                    // Calculate all forces of a node and the spring_forces
                    ImplementationComponents<Tag_, DataType_, WeightedKamadaKawai>::initializeForces(_spring_force_parameters, _coordinates, _spring_forces, node_forces);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements());
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(_spring_forces[i.index()]);
                        *i > result ? result = *i, _max_node = i.index() : 0;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > std::numeric_limits<DataType_>::epsilon()) 
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[_max_node].copy());
                        Scale<>::value(_spring_force_of_max_node, _step_width[_max_node] / result);
                        Sum<>::value(_coordinates[_max_node], _spring_force_of_max_node);
                    }

                    // Calculate statistic values
                    if ( (_statistic_step > 0) && (_number_of_iterations-1 == 0 || ((_number_of_iterations-1) % _statistic_step) == 0 ) )
                    {
                        DenseVector<DataType_> statistic_values(4, DataType_(0));
                        Statistics<Tag_>::value(_number_of_iterations-1, _weights_of_edges, result, _coordinates, statistic_values);
                        statistic_queue.push(statistic_values);
                    }

                    //std::cout << "_coordinates "<< _coordinates << std::endl;
                    _max_force = result;
                }

                DataType_ value(DataType_ & eps)
                {
                    // Increment number of iterations
                    _number_of_iterations++;

                    // previous_spring_forces
                    DenseMatrix<DataType_> _previous_spring_forces(_spring_forces.copy());

                    // previous_max_node
                    int _previous_max_node(_max_node);

                    // Calculate the new spring_forces, the maximal force, the node with maximal force and the new node_forces
                    DataType_ result(0);
                    ImplementationComponents<Tag_, DataType_, WeightedKamadaKawai>::calculateForces(_spring_forces, _max_node, _coordinates,
                    node_forces, _graph_distance[_previous_max_node], result);

                    // Calculate the new _step_width
                    DenseVectorRange<DataType_> previous_force(_previous_spring_forces[_max_node]);
                    DataType_ aux(Norm<vnt_l_two, true, Tag_>::value(_previous_spring_forces[_max_node]));
                    DenseVector<DataType_> current_force(_spring_forces[_max_node].copy());
                    for (typename DenseVector<DataType_>::ElementIterator i(previous_force.begin_elements()), i_end(previous_force.end_elements()), k(current_force.begin_elements());
                        i != i_end ; ++i, ++k)
                    {
                        *i /= aux;
                        *k /= result;
                    }
                    DataType_ prod( DotProduct<Tag_>::value(previous_force, current_force) );
                    if (prod > 0.8) _step_width[_max_node] *=_step_width_factor_1;
                    if ( (prod < -0.8) || (fabs(prod) < 0.2) ) _step_width[_max_node] *=_step_width_factor_2;

                    // Calculate the new positions by using resulting forces
                    DataType_ noise(1);
                    if (result > eps)
                    {
                        //_step_width[_max_node] > _noise_duration ? noise = (7.5 + ( rand() % 5 ) ) / 10 : noise = 1;
                        DataType_ delta(std::min(_step_width[_max_node], result));
                        delta *=noise;
                        Scale<>::value(current_force, delta);
                        Sum<>::value(_coordinates[_max_node], current_force);
                    }

                    // Calculate statistic values
                    if ( (_statistic_step > 0) && (_number_of_iterations-1 == 0 || ((_number_of_iterations-1) % _statistic_step) == 0 ) )
                    {
                        DenseVector<DataType_> statistic_values(4, DataType_(0));
                        Statistics<Tag_>::value(_number_of_iterations-1, _weights_of_edges, result, _coordinates, statistic_values);
                        statistic_queue.push(statistic_values);
                    }

                    _max_force = result;
                    //std::cout << "_coordinates in value "<< _coordinates << std::endl;

                    return result;
                }
        };
    }
}

#endif
