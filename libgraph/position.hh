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
#include <libgraph/node_distance_inverse.hh>
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
#include <libgraph/graph.hh>
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
                if (coordinates.columns() != neighbours.columns())
                   throw MatrixColumnsDoNotMatch(neighbours.columns(), coordinates.columns());

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

            Positions(DenseMatrix<DataType_> & coordinates, DenseVector<DataType_> & weights_of_nodes,
                    SparseMatrix<DataType_> & weights_of_edges) :
                _imp(new methods::Implementation<Tag_, DataType_, GraphTag_>(coordinates, weights_of_nodes, weights_of_edges))
            {
                if (coordinates.columns() != weights_of_edges.columns())
                    throw MatrixColumnsDoNotMatch(weights_of_edges.columns(), coordinates.columns());

                if (! weights_of_edges.square())
                    throw MatrixIsNotSquare(weights_of_edges.rows(), weights_of_edges.columns());

                if (weights_of_nodes.size() != coordinates.columns())
                    throw VectorSizeDoesNotMatch(weights_of_nodes.size(), coordinates.columns());
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
                    _step_width(DataType_(coordinates.columns() * edge_length / 20))
                {
                } 

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(NodeDistanceInverse<>::value(_coordinates));

                    // Having inv_square_dist, mask it and invert each element to get square_dist
                    SparseMatrix<DataType_> square_dist(_neighbours.columns(), _neighbours.rows());
                    for (typename Matrix<bool>::ConstElementIterator e(_neighbours.begin_non_zero_elements()),
                        e_end(_neighbours.end_non_zero_elements()); e != e_end ; ++e)
                    {
                        square_dist[e.row()][e.column()] = *e;
                    }

                    ElementProduct<>::value(square_dist, inv_square_dist);
                    ElementInverse<>::value(square_dist);

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum>::value(square_dist)));

                    // Calculate the same stuff for inv_square_dist
                    DenseVector<DataType_> * sum_vec_inv(new DenseVector<DataType_>(Reduction<rt_sum>::value(inv_square_dist)));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    Difference<>::value(diag, square_dist); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<>::value(DataType_(-1), square_dist);

                    // Calculating attractive_forces = (_coordinates * square_dist)
                    DenseMatrix<DataType_> attractive_forces(Product<>::value(_coordinates, square_dist));

                    // Calculating attractive_forces = (1/k^2 * attractive_forces)
                    Scale<>::value(DataType_(1 / _square_edge_length), attractive_forces);

                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<>::value(DataType_(-1), inv_square_dist);

                    // Calculating diff = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv->size(), *sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (_coordinates * diff)
                    DenseMatrix<DataType_> repulsive_forces(Product<>::value(_coordinates, diff));

                    // Calculating repulsive_forces = (k^2 * repulsive_forces)
                    Scale<>::value(_square_edge_length, repulsive_forces);

                    // Calculate the maximal force and the result forces
                    Sum<>::value(attractive_forces, repulsive_forces);
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true>::value(attractive_forces.column(i.index()));
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
                        result > eps ? *e = *e + _step_width / (resulting_forces[e.column()]) * *k : 0;
                    }

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
                DenseVector<DataType_> & _weights_of_nodes;

                ///  edge weight matrix - the weights of each edge
                SparseMatrix<DataType_> & _weights_of_edges;

                /// adjacence matrix - which nodes are neighbours?
                DenseMatrix<bool> _neighbours;

                /// parameter matrix for repulsive forces
                DenseMatrix<DataType_> _repulsive_force_parameter;

                /// parameter matrix for attractive forces
                DenseMatrix<DataType_> _attractive_force_parameter;

            public:
                friend class Positions<Tag_, DataType_, WeightedFruchtermanReingold>;

                Implementation(DenseMatrix<DataType_> & coordinates, DenseVector<DataType_> & weights_of_nodes,
                    SparseMatrix<DataType_> & weights_of_edges) :
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _neighbours(weights_of_edges.columns(), weights_of_edges.rows(), false),
                    _repulsive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), 0),
                    _attractive_force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), 0)
                {
                    // Calculate adjacence matrix
                    for (typename MutableMatrix<bool>::ElementIterator e(_neighbours.begin_elements()),
                        e_end(_neighbours.end_elements()); e != e_end ; ++e)
                    {
                        fabs(_weights_of_edges[e.row()][e.column()]) <= std::numeric_limits<DataType_>::epsilon() ?
                            *e = false : *e = true;
                    }

                    // Calculate parameter matrix for repulsive forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_repulsive_force_parameter.begin_elements()),
                        e_end(_repulsive_force_parameter.end_elements()); e != e_end ; ++e)
                    {
                        *e = (_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]) * (_weights_of_nodes[e.row()] * _weights_of_nodes[e.column()]);
                    }

                    // Calculate parameter matrix for attractive forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_attractive_force_parameter.begin_elements()),
                        e_end(_attractive_force_parameter.end_elements()), g(_weights_of_edges.begin_elements()); e != e_end ; ++e, ++g)
                    {
                        *e = *g * *g * *g * *g;
                    }
                }

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(NodeDistanceInverse<>::value(_coordinates));

                    // Having inv_square_dist, mask it and invert each element to get square_dist
                    DenseMatrix<DataType_> square_dist(*inv_square_dist.copy());
                    ElementProduct<>::value(square_dist, _neighbours);
                    ElementInverse<>::value(square_dist);

                    // Calculate inv_square_dist = Mul(inv_square_dist,  _repulsive_force_parameter)
                    ElementProduct<>::value(inv_square_dist, _repulsive_force_parameter);

                    // Calculate square_dist = Mul(square_dist,  _attractive_force_parameter)
                    ElementProduct<>::value(square_dist, _attractive_force_parameter);

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum>::value(square_dist)));

                    // Calculate the same stuff for inv_square_dist
                    DenseVector<DataType_> * sum_vec_inv(new DenseVector<DataType_>(Reduction<rt_sum>::value(inv_square_dist)));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    Difference<>::value(diag, square_dist);
                    Scale<>::value(DataType_(-1), square_dist);

                    // Calculating attractive_forces = (_coordinates * square_dist)
                    DenseMatrix<DataType_> attractive_forces(Product<>::value(_coordinates, square_dist));
                   
                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<>::value(DataType_(-1), inv_square_dist);

                    // Calculating inv_square_dist = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv->size(), *sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (_coordinates * inv_square_dist)
                    DenseMatrix<DataType_> repulsive_forces(Product<>::value(_coordinates, diff));

                    // Calculate the maximal force and the result forces
                    Sum<>::value(attractive_forces, repulsive_forces);
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true>::value(attractive_forces.column(i.index()));
                        result = std::max(result, *i);
                    }

                    // Calculate the new positions by using result forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_coordinates.begin_elements()),
                        e_end(_coordinates.end_elements()), k(attractive_forces.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        result > eps ? *e = *e + 1 / (10 * resulting_forces[e.column()]) * *k : 0;
                    }
                    return result;
                }

                void init()
                {
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
                    _step_widths(coordinates.columns(), DataType_(coordinates.columns() * edge_length / 20))

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
                    ElementProduct<>::value(_graph_distance, _graph_distance);
                    Scale<>::value(_square_edge_length, _graph_distance);

                    // Calculate 1 / _graph_distance
                    ElementInverse<>::value(_graph_distance);

                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<>::value(_coordinates));

                    // _spring_force_parameters = mul( _graph_distance, square_dist) -1
                    _spring_force_parameters = ElementProduct<>::value(*(_graph_distance.copy()), square_dist);
                    Sum<>::value(_spring_force_parameters, DataType_(-1));

                    // Calculate the single diagonal matrix containing the row sum vector of _spring_force_parameters
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum>::value(_spring_force_parameters)));

                    // Calculating square_dist = (_spring_force_parameters - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    square_dist = Difference<>::value(diag, *(_spring_force_parameters.copy())); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<>::value(DataType_(-1), square_dist);

                    // Calculating spring_forces = (_coordinates * square_dist)
                    _spring_forces = Product<>::value(_coordinates, square_dist);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements());
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true>::value(_spring_forces.column(i.index()));
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > std::numeric_limits<DataType_>::epsilon()) 
                    {
                        DenseVector<DataType_> _spring_force_of_max_node(_spring_forces.column(max_node).copy());
                        Scale<>::value(DataType_(1 / (result * 10)), _spring_force_of_max_node);
                        Sum<>::value(_coordinates.column(max_node), _spring_force_of_max_node);
                    }
                }


                DataType_ value(DataType_ & eps)
                {
                    // Force parameter vector of node with maximal force
                    DenseVector<DataType_> _spring_force_parameter_vector(_spring_force_parameters.column(max_node).copy());

                    // Calculate the difference between _previous_coordinates and the previous coordinates of max_node
                    DenseMatrix<DataType_> _previous_coordinates_difference(*(_previous_coordinates.copy()));
                    for (int i(0) ; i != _previous_coordinates_difference.columns() ; ++i)
                    {
                        Difference<Tag_>::value(_previous_coordinates_difference.column(i), _previous_coordinates.column(max_node));
                        Scale<>::value(DataType_(-1), _previous_coordinates_difference.column(i));
                    }

                    // Calculate the new _spring_force_parameters
                    DenseVector<DataType_> _spring_force_parameters_of_max_node(_spring_force_parameters.column(max_node));
                    for (int i(0) ; i != _spring_force_parameters_of_max_node.size() ; ++i)
                    {
                        DenseVector<DataType_> auxiliary(_coordinates.column(i).copy());
                        Difference<>::value(auxiliary, _coordinates.column(max_node));
                        DataType_ norm_of_auxiliary(Norm<vnt_l_two, false>::value(auxiliary));
                        norm_of_auxiliary == std::numeric_limits<DataType_>::epsilon() ? 
                        _spring_force_parameters_of_max_node[i] = DataType_(0) :
                        _spring_force_parameters_of_max_node[i] = norm_of_auxiliary;
                    }
                    ElementProduct<>::value(_spring_force_parameters_of_max_node, _graph_distance.column(max_node).copy());
                    Sum<>::value(DataType_(-1), _spring_force_parameters_of_max_node);
                    for (int i(0) ; i != _spring_force_parameters.columns() ; ++i)
                    {
                        _spring_force_parameters[max_node][i] = _spring_force_parameters[i][max_node];
                    }

                    // Calculate the difference between _coordinates and the coordinates of max_node
                    DenseMatrix<DataType_> _coordinates_difference(*(_coordinates.copy()));
                    for (int i(0) ; i != _coordinates_difference.columns() ; ++i)
                    {
                        Difference<>::value(_coordinates_difference.column(i), _coordinates.column(max_node));
                        Scale<>::value(DataType_(-1), _coordinates_difference.column(i));
                    }

                    // Calculate the force difference (subtract the previous force of max_node and add the current force of max_node)
                    for (int i(0) ; i != _coordinates_difference.columns() ; ++i)
                    {
                        if (i != max_node)
                        {
                            Scale<>::value(_spring_force_parameters[max_node][i], _coordinates_difference.column(i));
                            Scale<>::value(_spring_force_parameter_vector[i], _previous_coordinates_difference.column(i));
                            Sum<>::value(_coordinates_difference.column(max_node), _coordinates_difference.column(i));
                            Difference<>::value(_coordinates_difference.column(i), _previous_coordinates_difference.column(i));
                        }
                    }
                    Sum<>::value(_coordinates_difference.column(max_node), _spring_forces.column(max_node));
                    Scale<>::value(DataType_(-1), _coordinates_difference.column(max_node));

                    // previous_spring_force of max_node
                    DenseVector<DataType_> _previous_spring_force(_spring_forces.column(max_node).copy());

                    // Calculating _spring_forces = (_spring_forces + _force_difference)
                    Sum<>::value(_spring_forces, _coordinates_difference);

                    // previous_max_node
                    int _previous_max_node(max_node);

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true>::value(_spring_forces.column(i.index()));
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
                        for (int i(0); i < _spring_forces.rows(); i++)
                        {
                            if (fabs(_previous_spring_force[i] - _spring_forces[i][max_node] * (-1)) > fabs(0.3 * _previous_spring_force[i]))
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
                            DenseVector<DataType_> _spring_force_of_max_node(_spring_forces.column(max_node).copy());
                            Scale<>::value(DataType_(_step_widths[max_node] / result), _spring_force_of_max_node);
                            Sum<>::value(_coordinates.column(max_node), _spring_force_of_max_node);
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

                /// weight vector - the weights of each node
                DenseVector<DataType_> & _weights_of_nodes;

                /// edge weight matrix - the weights of each edge
                SparseMatrix<DataType_> & _weights_of_edges;

                /// parameter matrix for forces - square optimal distance
                DenseMatrix<DataType_> _force_parameter;

            public:
                friend class Positions<Tag_, DataType_, WeightedKamadaKawai>;
                
                Implementation(DenseMatrix<DataType_> & coordinates, DenseVector<DataType_> & weights_of_nodes,
                    SparseMatrix<DataType_> & weights_of_edges) :
                    _coordinates(coordinates),
                    _weights_of_nodes(weights_of_nodes),
                    _weights_of_edges(weights_of_edges),
                    _force_parameter(weights_of_edges.columns(), weights_of_edges.rows(), 0)
                {
                    // Calculate cost_matrix represented by _force_parameter to get previous nodes of unweigthed graph
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_force_parameter.begin_elements()),
                        e_end(_force_parameter.end_elements()); e != e_end ; ++e)
                    {
                        fabs(_weights_of_edges[e.row()][e.column()]) <= std::numeric_limits<DataType_>::epsilon() ?
                        *e = 2 * weights_of_edges.columns() : *e = 1;
                    }

                    // Initialize a matrix of previous_nodes
                    DenseMatrix<int> previous_nodes(weights_of_edges.columns(), weights_of_edges.rows(), 0);

                    // Calculate previous nodes
                    Dijkstra<DataType_>::value(_force_parameter, previous_nodes);

                    // Calculate a optimal distance, _force_parameter <- optimal distance
                    int row;
                    DataType_ costs;
                    for (typename MutableMatrix<int>::ElementIterator e(previous_nodes.begin_elements()),
                        e_end(previous_nodes.end_elements()); e != e_end ; ++e)
                    {
                        row = e.row();
                        costs = 0;
                        while (row != e.column())
                        {
                            costs = sqrt(_weights_of_nodes[row] * _weights_of_nodes[previous_nodes[row][e.column()]]) /
                            weights_of_edges[row][previous_nodes[row][e.column()]];
                            row = previous_nodes[row][e.column()];
                        }
                    _force_parameter[e.column()][e.row()] = costs;
                    }

                    ElementProduct<>::value (_force_parameter, _force_parameter);
                }

                void init()
                {
                    // TODO: Implement weighted KamadaKawai initialisaton
                }

                DataType_ value(DataType_ & eps)
                {
                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<>::value(_coordinates));

                    // dist = square_dist *2 
                    DenseMatrix<DataType_> dist(*square_dist.copy());
                    Scale<>::value(DataType_(2), dist);

                    // square_dist = square_dist + _force_parameter
                    Sum<>::value(square_dist, _force_parameter);

                    // square_dist = 1 / square_dist
                    ElementInverse<>::value(square_dist);

                    // Mul(square_dist, dist) -1
                    ElementProduct<>::value(square_dist, dist);
                    Sum<>::value(square_dist, DataType_(-1));

                    // Calculate the single diagonal matrix containing the row sum vector of square_dist
                    DenseVector<DataType_> * sum_vec(new DenseVector<DataType_>(Reduction<rt_sum>::value(square_dist)));

                    // Calculating square_dist = (square_dist - diag(sum_vec))
                    BandedMatrix<DataType_> diag(sum_vec->size(), *sum_vec);
                    Difference<>::value(diag, square_dist); /// \todo Switch to MD::value(Dense, const Banded)
                    Scale<>::value(DataType_(-1), square_dist);

                    // Calculating attractive_forces = (_coordinates * square_dist)
                    DenseMatrix<DataType_> attractive_forces(Product<>::value(_coordinates, square_dist));

                    // Calculate the resulting forces, the maximal force and the node with maximal force
                    DataType_ result(0);
                    int max_node(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0));
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = Norm<vnt_l_two, true>::value(attractive_forces.column(i.index()));
                        *i > result ? result = *i, max_node = i.index() : 0;
                    }

                    // Calculate the new positions by using resulting forces
                    if (result > eps)
                    {
                        for (int i(0) ; i != _coordinates.rows() ; ++i)
                        {
                            _coordinates[i][max_node] = _coordinates[i][max_node] +
                            1 / (10 * resulting_forces[max_node]) * attractive_forces[i][max_node];
                        }
                    }
                    return result;
                }
        };
    }
}

#endif
