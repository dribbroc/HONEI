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

#include <honei/libgraph/dijkstra.hh>
#include <honei/libgraph/node_distance.hh>
#include <honei/libgraph/position-impl.hh>
#include <honei/libgraph/breadth_first_search.hh>
#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/difference.hh>
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

    template <typename Tag_ , typename DataType_, typename GraphTag_> class Positions
    {
        private:
            /// Our implementation pattern.
            methods::Implementation<Tag_ ,DataType_, GraphTag_> * _imp;

        public:  
            Positions(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                _imp(new methods::Implementation<Tag_, DataType_, GraphTag_>(coordinates, weights_of_nodes, weights_of_edges))
            {
                if (coordinates.rows() != weights_of_edges.columns())
                    throw MatrixColumnsDoNotMatch(weights_of_edges.columns(), coordinates.columns());

                if (! weights_of_edges.square())
                    throw MatrixIsNotSquare(weights_of_edges.rows(), weights_of_edges.columns());

                if (weights_of_nodes.size() != coordinates.rows())
                    throw VectorSizeDoesNotMatch(weights_of_nodes.size(), coordinates.columns());
#if 0
                \todo Implement Trace
                if (Trace<>::value(weights_of_edges) > std::numeric_limits<DataType_>::epsilon())
                    throw MatrixHasNonVanishingTrace;
#endif
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

                /// previous max_Force
                DataType_ _previous_max_Force;

                /// ancestor of previous max_Force
                DataType_ _previous_max_Force_2;

                /// number of iterations
                unsigned long _number_of_iterations;

                /// the maximal force
                DataType_ _max_force;

                /// step width
                DataType_ _step_width;

                /// the _repulsive_force_range defines the range of the repulsive force
                DataType_ _repulsive_force_range;

            public:
                friend class Positions<Tag_, DataType_, WeightedFruchtermanReingold>;

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _repulsive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), DataType_(0)),
                    _attractive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), DataType_(0)),
                    _previous_max_Force(DataType_(0)),
                    _previous_max_Force_2(DataType_(0)),
                    _number_of_iterations(1),
                    _max_force(0),
                    _step_width(0),
                    _repulsive_force_range(DataType_(0))
                {
                // Calculate parameter matrix for repulsive forces and attractive forces, _repulsive_force_range and _step_width
                    typename MutableMatrix<DataType_>::ElementIterator f(_attractive_force_parameter.begin_elements());
                    typename Matrix<DataType_>::ConstElementIterator g(_weights_of_edges.begin_elements());
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_repulsive_force_parameter.begin_elements()),
                        e_end(_repulsive_force_parameter.end_elements()); e != e_end ; ++e, ++f, ++g)
                    {
                        DataType_ node_weight(_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]);
                        if (*g > std::numeric_limits<DataType_>::epsilon())
                        {
                            DataType_ length_of_edge(sqrt(node_weight) / *g);
                            if ((length_of_edge) > _repulsive_force_range) _repulsive_force_range = length_of_edge;
                        }
                        *e = node_weight * node_weight;
                        *f = *g * *g * *g * *g;
                    }
                    _repulsive_force_range *= _weights_of_edges.rows();
                    _step_width = _repulsive_force_range / 20;
                }

                Implementation(AbstractGraph<DataType_> & graph, DataType_ edgeLength) :
                    _coordinates(*graph.coordinates()),
                    _weights_of_nodes(*graph.nodeWeights()),
                    _weights_of_edges(*graph.edges()),
                    _repulsive_force_parameter(_weights_of_edges.columns(), _weights_of_edges.rows(), DataType_(0)),
                    _attractive_force_parameter(_weights_of_edges.columns(), _weights_of_edges.rows(), DataType_(0)),
                    _previous_max_Force(DataType_(0)),
                    _previous_max_Force_2(DataType_(0)),
                    _number_of_iterations(1),
                    _max_force(0),
                    _step_width(0),
                    _repulsive_force_range(DataType_(0))
                {
                    // Calculate parameter matrix for repulsive forces and attractive forces, _repulsive_force_range and _step_width
                    typename MutableMatrix<DataType_>::ElementIterator f(_attractive_force_parameter.begin_elements());
                    typename Matrix<DataType_>::ConstElementIterator g(_weights_of_edges.begin_elements());
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_repulsive_force_parameter.begin_elements()),
                        e_end(_repulsive_force_parameter.end_elements()); e != e_end ; ++e, ++f, ++g)
                    {
                        DataType_ node_weight(_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]);
                        if (*g > std::numeric_limits<DataType_>::epsilon())
                        {
                            DataType_ length_of_edge(sqrt(node_weight) / *g);
                            if ((length_of_edge) > _repulsive_force_range) _repulsive_force_range = length_of_edge;
                        }
                        *e = (graph.sameTimeslice(e.row(), e.column())) ?
                           node_weight * node_weight: 0;
                        *f = *g * *g * *g * *g;
                    }
                    _repulsive_force_range *= _weights_of_edges.rows();
                    _step_width = _repulsive_force_range / 20;
                }

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2) and _square_dist = (d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0));
                    DenseMatrix<DataType_> square_dist(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0));
                    NodeDistance<Tag_>::value(_coordinates, _weights_of_edges, square_dist, inv_square_dist, _repulsive_force_range);

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
                    if ((_number_of_iterations > 2) && (fabs(result - _previous_max_Force_2) <= fabs(0.2 * _previous_max_Force_2)))
                    {
                        if (_step_width > _repulsive_force_range / (_weights_of_edges.rows() * 20)) _step_width *= 0.5;
                        if (_step_width <= _repulsive_force_range  / (_weights_of_edges.rows() * 20)) _step_width *= 0.995;
                    }

                    // Calculate the new previous forces
                    _previous_max_Force_2 = _previous_max_Force;
                    _previous_max_Force = result;

                    // Calculate the new positions by using result forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_coordinates.begin_elements()),
                        e_end(_coordinates.end_elements()), k(attractive_forces.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        result > eps ? *e = *e + _step_width / (resulting_forces[e.row()]) * *k : 0;
                    }

                    _number_of_iterations++;
                    _max_force = result;

                    return result;
                }

                void init()
                {

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

                /// vector of step widths
                DenseVector<DataType_> _step_widths;

            public:
                friend class Positions<Tag_, DataType_, WeightedKamadaKawai>;

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
                    _step_widths(coordinates.rows(), DataType_(0))
                {
                    // Using BFS to calculate the graph distance matrix
                    if (! BreadthFirstSearch<>::value(_graph_distance, _weights_of_nodes, _weights_of_edges))
                        throw GraphError("Graph has to be coherently");
                }

                Implementation(AbstractGraph<DataType_> & graph, DataType_ edge_length) :
                    _coordinates(*graph.coordinates()),
                    _weights_of_nodes(*graph.nodeWeights()),
                    _weights_of_edges(*graph.edges()),
                    _graph_distance(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0)),
                    _spring_forces(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _max_force(0),
                    _number_of_iterations(1),
                    _max_node(0),
                    node_forces(_coordinates.rows() * _coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _step_widths(_coordinates.rows(), DataType_(0))
                {
                    // Using BFS to calculate the graph distance matrix
                    if (! BreadthFirstSearch<>::value(_graph_distance, _weights_of_nodes, _weights_of_edges, graph))
                        throw GraphError("Graph has to be coherently");
                }

                void init()
                {
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

                    // Fill vector of _step_widths with stepwidth
                    for (typename DenseVector<DataType_>::ElementIterator i(_step_widths.begin_elements()), i_end(_step_widths.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = stepwidth;
                    }

                    // Calculate square _graph_distance
                    ElementProduct<Tag_>::value(_graph_distance, _graph_distance);

                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<Tag_>::value(_coordinates));

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
                        Scale<>::value(_spring_force_of_max_node, _step_widths[_max_node] / result);
                        Sum<>::value(_coordinates[_max_node], _spring_force_of_max_node);
                    }

                    //std::cout << "_coordinates "<< _coordinates << std::endl;
                    _max_force = result;

                }

                DataType_ value(DataType_ & eps)
                {
                    // Increment number of iterations
                    _number_of_iterations++;

                    // previous_spring_force of _max_node
                    DenseVector<DataType_> _previous_spring_force(_spring_forces[_max_node].copy());
                    Scale<>::value(_previous_spring_force, DataType_(1 / Norm<vnt_l_two, true, Tag_>::value(_previous_spring_force)));

                    // previous_max_node
                    int _previous_max_node(_max_node);

                    // Calculate the new spring_forces, the maximal force, the node with maximal force and the new node_forces
                    DataType_ result(0);
                    ImplementationComponents<Tag_, DataType_, WeightedKamadaKawai>::calculateForces(_spring_forces, _max_node, _coordinates,
                    node_forces, _graph_distance[_previous_max_node], result);

                    // Calculate the new step_width
                    if (_max_node == _previous_max_node)
                    {
                        bool reducing_condition(true);
                        for (int i(0); i < _spring_forces.columns(); i++)
                        {
                            if (fabs(_previous_spring_force[i] - _spring_forces[_max_node][i] * (-1 / result)) > fabs(0.3 * _previous_spring_force[i]))
                            {
                                reducing_condition = false;
                                break;
                            }
                        }
                        if (reducing_condition) _step_widths[_max_node] = _step_widths[_max_node] / 2;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > eps)
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[_max_node].copy());
                        Scale<>::value(_spring_force_of_max_node, DataType_(_step_widths[_max_node] / result));
                        Sum<>::value(_coordinates[_max_node], _spring_force_of_max_node);
                    }

                    _max_force = result;
                    //std::cout << "_coordinates in value "<< _coordinates << std::endl;

                    return result;
                }
        };
    }
}

#endif
