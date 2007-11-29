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

#include <libgraph/dijkstra.hh>
#include <libgraph/node_distance.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/difference.hh>
#include <libla/element_inverse.hh>
#include <libla/element_product.hh>
#include <libla/matrix_error.hh>
#include <libla/product.hh>
#include <libla/reduction.hh>
#include <libla/sum.hh>
#include <libla/scale.hh>
#include <libla/vector.hh>
#include <libgraph/abstract_graph.hh>
#include <libgraph/graph_error.hh>
#include <libutil/tags.hh>

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
         * Implementation tag for Fruchterman-Reingold method.
         */
        struct FruchtermanReingold
        {
        };

        /**
         * Implementation tag for Kamada-Kawai method.
         */
        struct KamadaKawai
        {
        };

        /**
         * Implementation tag for weighted Fruchterman-Reingold method.
         */
        struct WeightedFruchtermanReingold
        {
        };

        /**
         * Implementation tag for weighted Kamada-Kawai method.
         */
        struct WeightedKamadaKawai
        {
        };

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

            Positions(DenseMatrix<DataType_> & coordinates, const SparseMatrix<bool> & neighbours, DataType_ edge_length) :
                _imp(new methods::Implementation<Tag_, DataType_, GraphTag_>(coordinates, neighbours, edge_length))
            {
                if (coordinates.rows() != neighbours.columns())
                   throw MatrixColumnsDoNotMatch(neighbours.columns(), coordinates.rows());

                if (! neighbours.square())
                    throw MatrixIsNotSquare(neighbours.rows(), neighbours.columns());
#if 0
                \todo Implement Trace
                if (Trace<>::value(neighbours) > false)
                    throw MatrixHasNonVanishingTrace;
#endif
                if (edge_length <= std::numeric_limits<DataType_>::epsilon())
                    throw GraphError("Edge length must be positive");
            }

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

            void update(DataType_ eps, int timeout)
            {
                _imp->init();
                do
                {
                    --timeout;
                }
                while ((timeout > 0) && (_imp->value(eps) > eps));
            }

            DenseMatrix<DataType_> & coordinates()
            {
                return _imp->_coordinates;
            }
    };

    namespace methods
    {
        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, FruchtermanReingold>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                /// adjacence matrix - which nodes are neighbours?
                const SparseMatrix<bool> & _neighbours;

                /// a scalar factor which determinates how much influence the distance has while calculating the forces.
                const DataType_ _edge_length;

                /// a scalar factor which determinates how much influence the distance has while calculating the forces.
                const DataType_ _square_edge_length;

                /// previous max_Force
                DataType_ _previous_max_Force;

                /// ancestor of previous max_Force
                DataType_ _previous_max_Force_2;

                /// number of iterations
                int _number_of_iterations;

                /// step width
                DataType_ _step_width;

                /// the _repulsive_force_range defines the range of the repulsive force
                DataType_ _repulsive_force_range;

            public:
                friend class Positions<Tag_, DataType_, FruchtermanReingold>;

                Implementation(DenseMatrix<DataType_> & coordinates, const SparseMatrix<bool> & neighbours, const DataType_ edge_length) :
                    _coordinates(coordinates),
                    _neighbours(neighbours),
                    _edge_length(edge_length),
                    _square_edge_length(edge_length * edge_length),
                    _previous_max_Force(0),
                    _previous_max_Force_2(0),
                    _number_of_iterations(1),
                    _step_width(DataType_(coordinates.rows() * edge_length / 20)),
                    _repulsive_force_range(coordinates.rows() * edge_length)
                {
                } 

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2) and _square_dist = (d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(_neighbours.columns(), _neighbours.rows(), 0);
                    SparseMatrix<DataType_> square_dist(_neighbours.columns(), _neighbours.rows());
                    NodeDistance<Tag_>::value(_coordinates, _neighbours, square_dist, inv_square_dist, _repulsive_force_range);

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(square_dist)));

                    // Calculate the same stuff for inv_square_dist
                    DenseVector<DataType_> * sum_vec_inv(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(inv_square_dist)));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    Difference<Tag_>::value(diag, square_dist); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<Tag_>::value(DataType_(-1), square_dist);

                    // Calculating attractive_forces = (square_dist * _coordinates)
                    DenseMatrix<DataType_> attractive_forces(Product<Tag_>::value(square_dist, _coordinates));

                    // Calculating attractive_forces = (1/k^2 * attractive_forces)
                    Scale<Tag_>::value(DataType_(1 / _square_edge_length), attractive_forces);

                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<Tag_>::value(DataType_(-1), inv_square_dist);

                    // Calculating diff = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv->size(), *sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<Tag_>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (diff * _coordinates)
                    DenseMatrix<DataType_> repulsive_forces(Product<Tag_>::value(diff, _coordinates));

                    // Calculating repulsive_forces = (k^2 * repulsive_forces)
                    Scale<Tag_>::value(_square_edge_length, repulsive_forces);

                    // Calculate the maximal force and the result forces
                    Sum<Tag_>::value(attractive_forces, repulsive_forces);
                    DataType_ result(0);
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
                        if (_step_width > _edge_length / 20) _step_width *= 0.5;
                        if (_step_width <= _edge_length / 20) _step_width *= 0.995;
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

                    return result;
                }

                void init()
                {
                }

        };

        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, WeightedFruchtermanReingold>
        {
            private:
                ///  position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                ///  weight vector - the weights of each node
                const DenseVector<DataType_> & _weights_of_nodes;

                ///  edge weight matrix - the weights of each edge
                const SparseMatrix<DataType_> & _weights_of_edges;

                /// parameter matrix for repulsive forces
                DenseMatrix<DataType_> _repulsive_force_parameter;

                /// parameter matrix for attractive forces
                DenseMatrix<DataType_> _attractive_force_parameter;

                /// previous max_Force
                DataType_ _previous_max_Force;

                /// ancestor of previous max_Force
                DataType_ _previous_max_Force_2;

                /// number of iterations
                int _number_of_iterations;

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
                    _repulsive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), 0),
                    _attractive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), 0),
                    _previous_max_Force(0),
                    _previous_max_Force_2(0),
                    _number_of_iterations(1),
                    _step_width(0),
                    _repulsive_force_range(0)
                {
                }

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2) and _square_dist = (d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(_weights_of_edges.rows(), _weights_of_edges.columns(), 0);
                    SparseMatrix<DataType_> square_dist(_weights_of_edges.rows(), _weights_of_edges.columns());
                    NodeDistance<Tag_>::value(_coordinates, _weights_of_edges, square_dist, inv_square_dist, _repulsive_force_range);

                    // Calculate inv_square_dist = Mul(inv_square_dist,  _repulsive_force_parameter)
                    ElementProduct<Tag_>::value(inv_square_dist, _repulsive_force_parameter);

                    // Calculate square_dist = Mul(square_dist,  _attractive_force_parameter)
                    ElementProduct<Tag_>::value(square_dist, _attractive_force_parameter);

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(square_dist)));

                    // Calculate the same stuff for inv_square_dist
                    DenseVector<DataType_> * sum_vec_inv(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(inv_square_dist)));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    Difference<Tag_>::value(diag, square_dist);
                    Scale<Tag_>::value(DataType_(-1), square_dist);

                    // Calculating attractive_forces = (square_dist * _coordinates)
                    DenseMatrix<DataType_> attractive_forces(Product<Tag_>::value(square_dist, _coordinates));
                   
                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<Tag_>::value(DataType_(-1), inv_square_dist);

                    // Calculating diff = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv->size(), *sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<Tag_>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (diff * _coordinates)
                    DenseMatrix<DataType_> repulsive_forces(Product<Tag_>::value(diff, _coordinates));

                    // Calculate the maximal force and the result forces
                    Sum<Tag_>::value(attractive_forces, repulsive_forces);
                    DataType_ result(0);
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
                        if (_step_width > _repulsive_force_range / 20) _step_width *= 0.5;
                        if (_step_width <= _repulsive_force_range  / 20) _step_width *= 0.995;
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

                    return result;
                }

                void init()
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

        };

        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, KamadaKawai>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                /// previous position matrix
                DenseMatrix<DataType_> & _previous_coordinates;

                /// adjacence matrix - which nodes are neighbours?
                const SparseMatrix<bool> & _neighbours;

                /// a scalar factor which determinates how much influence the distance has while calculating the forces.
                const DataType_ _square_edge_length;

                /// graph distance matrix
                DenseMatrix<DataType_> _graph_distance;

                /// spring forces
                DenseMatrix<DataType_> _spring_forces;

                /// parameter matrix for spring forces
                DenseMatrix<DataType_> _spring_force_parameters;

                /// index of node with maximal force
                int max_node ;

                /// vector of step widths
                DenseVector<DataType_> _step_widths;

            public:
                friend class Positions<Tag_, DataType_, KamadaKawai>;

                Implementation(DenseMatrix<DataType_> & coordinates, const SparseMatrix<bool> & neighbours, const DataType_ edge_length) :
                    _previous_coordinates(coordinates),
                    _coordinates(*(coordinates.copy())),
                    _neighbours(neighbours),
                    _square_edge_length(DataType_(edge_length * edge_length)),
                    _graph_distance(neighbours.columns(), neighbours.rows(), DataType_(2) * neighbours.columns()),
                    _spring_force_parameters(neighbours.columns(), neighbours.rows(), DataType_(0)),
                    _spring_forces(coordinates.columns(), coordinates.rows(), DataType_(0)),
                    max_node(0),
                    _step_widths(coordinates.rows(), DataType_(coordinates.rows() * edge_length / 20))

                {

                }

                void init()
                {
                    // Use adjacency matrix to calculate the _graph_distance
                    for (typename Matrix<bool>::ConstElementIterator e(_neighbours.begin_non_zero_elements()),
                        e_end(_neighbours.end_non_zero_elements()) ; e != e_end ; ++e)
                    {
                        _graph_distance[e.row()][e.column()] = 1;
                    }

                    // Using Dijkstra to calculate the graph distance matrix
                    Dijkstra<DataType_>::value(_graph_distance);

                    // Mul(_graph_distance, _graph_distance) * _edge_length^2
                    ElementProduct<Tag_>::value(_graph_distance, _graph_distance);
                    Scale<Tag_>::value(_square_edge_length, _graph_distance);

                    // Calculate 1 / _graph_distance
                    ElementInverse<Tag_>::value(_graph_distance);

                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<Tag_>::value(_coordinates));

                    // _spring_force_parameters = mul( _graph_distance, square_dist) -1
                    _spring_force_parameters = ElementProduct<Tag_>::value(*(_graph_distance.copy()), square_dist);
                    Sum<Tag_>::value(_spring_force_parameters, DataType_(-1));


                    // Calculate the single diagonal matrix containing the row sum vector of _spring_force_parameters
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(_spring_force_parameters)));

                    // Calculating square_dist = (_spring_force_parameters - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    square_dist = Difference<Tag_>::value(diag, *(_spring_force_parameters.copy())); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<Tag_>::value(DataType_(-1), square_dist);

                    // Calculating spring_forces = (_coordinates * square_dist)
                    _spring_forces = Product<Tag_>::value(square_dist, _coordinates);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements());
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(_spring_forces[i.index()]);
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > std::numeric_limits<DataType_>::epsilon()) 
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[max_node].copy());
                        Scale<Tag_>::value(_step_widths[max_node] / result, _spring_force_of_max_node);
                        Sum<Tag_>::value(_coordinates[max_node], _spring_force_of_max_node);
                    }
                }


                DataType_ value(DataType_ & eps)
                {
                    // Force parameter vector of node with maximal force
                    DenseVector<DataType_> _spring_force_parameter_vector(_spring_force_parameters[max_node].copy());

                    // Calculate the difference between _previous_coordinates and the previous coordinates of max_node
                    DenseMatrix<DataType_> _previous_coordinates_difference(*(_previous_coordinates.copy()));
                    for (int i(0) ; i != _previous_coordinates_difference.rows() ; ++i)
                    {
                        Difference<Tag_>::value(_previous_coordinates_difference[i], _previous_coordinates[max_node]);
                        Scale<Tag_>::value(DataType_(-1), _previous_coordinates_difference[i]);
                    }

                    // Calculate the new _spring_force_parameters
                    DenseVectorRange<DataType_> _spring_force_parameters_of_max_node(_spring_force_parameters[max_node]);
                    for (int i(0) ; i != _spring_force_parameters_of_max_node.size() ; ++i)
                    {
                        DenseVector<DataType_> auxiliary(_coordinates[i].copy());
                        Difference<Tag_>::value(auxiliary, _coordinates[max_node]);
                        DataType_ norm_of_auxiliary(Norm<vnt_l_two, false, Tag_>::value(auxiliary));
                        norm_of_auxiliary > std::numeric_limits<DataType_>::epsilon() ? 
                        _spring_force_parameters_of_max_node[i] = norm_of_auxiliary:
                        _spring_force_parameters_of_max_node[i] = DataType_(0);
                    }

                    ElementProduct<Tag_>::value(_spring_force_parameters_of_max_node, _graph_distance[max_node].copy());
                    Sum<Tag_>::value(DataType_(-1), _spring_force_parameters_of_max_node);
                    for (int i(0) ; i != _spring_force_parameters.rows() ; ++i)
                    {
                        _spring_force_parameters[i][max_node] = _spring_force_parameters[max_node][i];
                    }

                    // Calculate the difference between _coordinates and the coordinates of max_node
                    DenseMatrix<DataType_> _coordinates_difference(*(_coordinates.copy()));
                    for (int i(0) ; i != _coordinates_difference.rows() ; ++i)
                    {
                        Difference<Tag_>::value(_coordinates_difference[i], _coordinates[max_node]);
                        Scale<Tag_>::value(DataType_(-1), _coordinates_difference[i]);
                    }

                    // Calculate the force difference (subtract the previous force of max_node and add the current force of max_node)
                    for (int i(0) ; i != _coordinates_difference.rows() ; ++i)
                    {
                        if (i != max_node)
                        {
                            Scale<Tag_>::value(_spring_force_parameters[max_node][i], _coordinates_difference[i]);
                            Scale<Tag_>::value(_spring_force_parameter_vector[i], _previous_coordinates_difference[i]);
                            Sum<Tag_>::value(_coordinates_difference[max_node], _coordinates_difference[i]);
                            Difference<Tag_>::value(_coordinates_difference[i], _previous_coordinates_difference[i]);
                        }
                    }
                    Sum<Tag_>::value(_coordinates_difference[max_node], _spring_forces[max_node]);
                    Scale<Tag_>::value(DataType_(-1), _coordinates_difference[max_node]);

                    // previous_spring_force of max_node
                    DenseVector<DataType_> _previous_spring_force(_spring_forces[max_node].copy());
                    Scale<Tag_>::value(DataType_(1 / Norm<vnt_l_two, true, Tag_>::value(_previous_spring_force)), _previous_spring_force);

                    // Calculating _spring_forces = (_spring_forces + _force_difference)
                    Sum<Tag_>::value(_spring_forces, _coordinates_difference);

                    // previous_max_node
                    int _previous_max_node(max_node);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements());
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(_spring_forces[i.index()]);
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // _previous_coordinates = _coordinates
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_previous_coordinates.begin_elements()),
                        e_end(_previous_coordinates.end_elements()), k(_coordinates.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        *e = *k;
                    }

                    // Calculate the new step_width
                    if (max_node == _previous_max_node)
                    {
                        bool reducing_condition(true);
                        for (int i(0); i < _spring_forces.columns(); i++)
                        {
                            if (fabs(_previous_spring_force[i] - _spring_forces[max_node][i] * (-1 / result)) > fabs(0.3 * _previous_spring_force[i]))
                            {
                                reducing_condition = false;
                                break;
                            }
                        }

                        if (reducing_condition) _step_widths[max_node] = _step_widths[max_node] / 2;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > eps)
                        {
                            DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[max_node].copy());
                            Scale<Tag_>::value(DataType_(_step_widths[max_node] / result), _spring_force_of_max_node);
                            Sum<Tag_>::value(_coordinates[max_node], _spring_force_of_max_node);
                        }
                    return result;
                }
        };

        template <typename Tag_, typename DataType_>
        class Implementation<Tag_, DataType_, WeightedKamadaKawai>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                /// previous position matrix
                DenseMatrix<DataType_> & _previous_coordinates;

                /// weight vector - the weights of each node
                const DenseVector<DataType_> & _weights_of_nodes;

                /// edge weight matrix - the weights of each edge
                const SparseMatrix<DataType_> & _weights_of_edges;

                /// graph distance matrix
                DenseMatrix<DataType_> _graph_distance;

                /// spring forces
                DenseMatrix<DataType_> _spring_forces;

                /// parameter matrix for spring forces
                DenseMatrix<DataType_> _spring_force_parameters;

                /// index of node with maximal force
                int max_node;

                /// vector of step widths
                DenseVector<DataType_> _step_widths;

                /// Reference to graph object, if this method is used with graph types. 
                AbstractGraph<DataType_> * _graph; 
                
            public:
                friend class Positions<Tag_, DataType_, WeightedKamadaKawai>;
                
                Implementation(DenseMatrix<DataType_> & coordinates, const DenseVector<DataType_> & weights_of_nodes,
                    const SparseMatrix<DataType_> & weights_of_edges) :
                    _previous_coordinates(*(coordinates.copy())),
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _graph_distance(weights_of_edges.rows(), weights_of_edges.columns(), DataType_(0)),
                    _spring_forces(coordinates.rows(), coordinates.columns(), DataType_(0)),
                    _spring_force_parameters(weights_of_edges.rows(), weights_of_edges.rows(), DataType_(0)),
                    max_node(0),
                    _step_widths(coordinates.rows(), DataType_(0)),
                    _graph(0)
                {
                }

                Implementation(AbstractGraph<DataType_> & graph, int edge_length) :
                    _coordinates(*graph.coordinates()),
                    _previous_coordinates(*(graph.coordinates()->copy())),
                    _weights_of_nodes(*graph.nodeWeights()),
                    _weights_of_edges(*graph.edges()),
                    _graph_distance(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0)),
                    _spring_forces(_coordinates.rows(), _coordinates.columns(), DataType_(0)),
                    _spring_force_parameters(_weights_of_edges.rows(), _weights_of_edges.columns(), DataType_(0)),
                    max_node(0),
                    _step_widths(_coordinates.rows(), DataType_(0)),
                    _graph(&graph)
                {
                }

                void init()
                {
                    // Calculate cost_matrix represented by _graph_distance to get previous nodes of unweigthed graph
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_graph_distance.begin_elements()),
                        e_end(_graph_distance.end_elements()); e != e_end ; ++e)
                    {
                        fabs(_weights_of_edges[e.row()][e.column()]) <= std::numeric_limits<DataType_>::epsilon() ?
                        *e = 2 * _weights_of_edges.columns() : *e = 1;
                    }

                    // Initialize a matrix of previous_nodes
                    DenseMatrix<int> previous_nodes(_weights_of_edges.columns(), _weights_of_edges.rows(), 0);

                    // Calculate previous nodes
                    Dijkstra<DataType_>::value(_graph_distance, previous_nodes);

                    // Calculate a optimal distance
                    int row;
                    DataType_ costs;
                    for (typename MutableMatrix<int>::ElementIterator e(previous_nodes.begin_elements()),
                        e_end(previous_nodes.end_elements()); e != e_end ; ++e)
                    {
                        row = e.row();
                        costs = 0;
                        while (row != e.column())
                        {
                            costs += (_graph == 0) || _graph->sameTimeslice(row, e.column()) ?
                            sqrt(_weights_of_nodes[row] * _weights_of_nodes[previous_nodes[row][e.column()]]) /
                            _weights_of_edges[row][previous_nodes[row][e.column()]] :
                            0;
                            row = previous_nodes[row][e.column()];
                        }
                    _graph_distance[e.column()][e.row()] = costs;
                    }

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
                    _spring_force_parameters = Sum<Tag_>::value(*(_graph_distance.copy()), square_dist);
                    ElementInverse<Tag_>::value(_spring_force_parameters);

                    // Calculate square_dist = square_dist * 2
                    Scale<Tag_>::value(DataType_(2), square_dist);

                    // _spring_force_parameters = _spring_force_parameters * square_dist -1
                    ElementProduct<Tag_>::value(_spring_force_parameters, square_dist);
                    Sum<Tag_>::value(_spring_force_parameters, DataType_(-1));


                    // Calculate the single diagonal matrix containing the row sum vector of _spring_force_parameters
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum, Tag_>::value(_spring_force_parameters)));

                    // Calculating square_dist = (_spring_force_parameters - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    square_dist = Difference<Tag_>::value(diag, *(_spring_force_parameters.copy())); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<Tag_>::value(DataType_(-1), square_dist);

                    // Calculating spring_forces = (_coordinates * square_dist)
                    _spring_forces = Product<Tag_>::value(square_dist, _coordinates);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements());
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(_spring_forces[i.index()]);
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > std::numeric_limits<DataType_>::epsilon()) 
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[max_node].copy());
                        Scale<Tag_>::value(_step_widths[max_node] / result, _spring_force_of_max_node);
                        Sum<Tag_>::value(_coordinates[max_node], _spring_force_of_max_node);
                    }

                }

                DataType_ value(DataType_ & eps)
                {
                    // Force parameter vector of node with maximal force
                    DenseVector<DataType_> _spring_force_parameter_vector(_spring_force_parameters[max_node].copy());

                    // Calculate the difference between _previous_coordinates and the previous coordinates of max_node
                    DenseMatrix<DataType_> _previous_coordinates_difference(*(_previous_coordinates.copy()));
                    for (int i(0) ; i != _previous_coordinates_difference.rows() ; ++i)
                    {
                        Difference<Tag_>::value(_previous_coordinates_difference[i], _previous_coordinates[max_node]);
                        Scale<Tag_>::value(DataType_(-1), _previous_coordinates_difference[i]);
                    }

                    // Calculate the new _spring_force_parameters
                    DenseVectorRange<DataType_> _spring_force_parameters_of_max_node(_spring_force_parameters[max_node]);
                    for (int i(0) ; i != _spring_force_parameters_of_max_node.size() ; ++i)
                    {
                        DenseVector<DataType_> auxiliary(_coordinates[i].copy());
                        Difference<Tag_>::value(auxiliary, _coordinates[max_node]);
                        DataType_ norm_of_auxiliary(Norm<vnt_l_two, false, Tag_>::value(auxiliary));
                        norm_of_auxiliary > std::numeric_limits<DataType_>::epsilon() ? 
                        _spring_force_parameters_of_max_node[i] = norm_of_auxiliary:
                        _spring_force_parameters_of_max_node[i] = DataType_(0);
                    }

                    DenseVector<DataType_> auxiliary2(_spring_force_parameters_of_max_node.copy());
                    Sum<Tag_>::value(auxiliary2, _graph_distance[max_node].copy());
                    ElementInverse<Tag_>::value(auxiliary2);
                    Scale<Tag_>::value(DataType_(2), _spring_force_parameters_of_max_node);
                    ElementProduct<Tag_>::value(_spring_force_parameters_of_max_node, auxiliary2);
                    Sum<Tag_>::value(DataType_(-1), _spring_force_parameters_of_max_node);
                    for (int i(0) ; i != _spring_force_parameters.rows() ; ++i)
                    {
                        _spring_force_parameters[i][max_node] = _spring_force_parameters[max_node][i];
                    }

                    // Calculate the difference between _coordinates and the coordinates of max_node
                    DenseMatrix<DataType_> _coordinates_difference(*(_coordinates.copy()));
                    for (int i(0) ; i != _coordinates_difference.rows() ; ++i)
                    {
                        Difference<Tag_>::value(_coordinates_difference[i], _coordinates[max_node]);
                        Scale<Tag_>::value(DataType_(-1), _coordinates_difference[i]);
                    }

                    // Calculate the force difference (subtract the previous force of max_node and add the current force of max_node)
                    for (int i(0) ; i != _coordinates_difference.rows() ; ++i)
                    {
                        if (i != max_node)
                        {
                            Scale<Tag_>::value(_spring_force_parameters[max_node][i], _coordinates_difference[i]);
                            Scale<Tag_>::value(_spring_force_parameter_vector[i], _previous_coordinates_difference[i]);
                            Sum<Tag_>::value(_coordinates_difference[max_node], _coordinates_difference[i]);
                            Difference<Tag_>::value(_coordinates_difference[i], _previous_coordinates_difference[i]);
                        }
                    }
                    Sum<Tag_>::value(_coordinates_difference[max_node], _spring_forces[max_node]);
                    Scale<Tag_>::value(DataType_(-1), _coordinates_difference[max_node]);

                    // previous_spring_force of max_node
                    DenseVector<DataType_> _previous_spring_force(_spring_forces[max_node].copy());
                    Scale<Tag_>::value(DataType_(1 / Norm<vnt_l_two, true, Tag_>::value(_previous_spring_force)), _previous_spring_force);

                    // Calculating _spring_forces = (_spring_forces + _force_difference)
                    Sum<Tag_>::value(_spring_forces, _coordinates_difference);

                    // previous_max_node
                    int _previous_max_node(max_node);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.rows(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true, Tag_>::value(_spring_forces[i.index()]);
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // _previous_coordinates = _coordinates
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_previous_coordinates.begin_elements()),
                        e_end(_previous_coordinates.end_elements()), k(_coordinates.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        *e = *k;
                    }

                    // Calculate the new step_width
                    if ((max_node == _previous_max_node) && (result > std::numeric_limits<DataType_>::epsilon()))
                    {
                        bool reducing_condition(true);
                        for (int i(0); i < _spring_forces.columns(); i++)
                        {
                            if (fabs(_previous_spring_force[i] - _spring_forces[max_node][i] * (-1 / result)) > fabs(0.3 * _previous_spring_force[i]))
                            {
                                reducing_condition = false;
                                break;
                            }
                        }

                        if (reducing_condition) _step_widths[max_node] = _step_widths[max_node] / 2;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > eps)
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces[max_node].copy());
                        Scale<Tag_>::value(DataType_(_step_widths[max_node] / result), _spring_force_of_max_node);
                        Sum<Tag_>::value(_coordinates[max_node], _spring_force_of_max_node);
                    }
                    return result;
                }
        };
    }
}

#endif
