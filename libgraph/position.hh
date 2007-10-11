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
        template <typename DataType_, typename Tag_>
        class Implementation;

    }

    template <typename DataType_, typename GraphTag_> class Positions
    {
        private:
            /// Our implementation pattern.
            methods::Implementation<DataType_, GraphTag_> * _imp;

        public:

            Positions(DenseMatrix<DataType_> & coordinates, const DenseMatrix<bool> & neighbours, DataType_ edge_length) :
                _imp(new methods::Implementation<DataType_, GraphTag_>(coordinates, neighbours, edge_length))
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
                _imp(new methods::Implementation<DataType_, GraphTag_>(coordinates, weights_of_nodes, weights_of_edges))
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
        template <typename DataType_>
        class Implementation<DataType_, FruchtermanReingold>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                /// adjacence matrix - which nodes are neighbours?
                const DenseMatrix<bool> & _neighbours;

                /// a scalar factor which determinates how much influence the distance has while calculating the forces.
                DataType_ _edge_length;

            public:
                friend class Positions<DataType_, FruchtermanReingold>;

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseMatrix<bool> & neighbours, DataType_ edge_length) :
                    _coordinates(coordinates),
                    _neighbours(neighbours),
                    _edge_length(edge_length)
                {
                } 

                DataType_ value(const DataType_ & eps)
                {
                    // Calculate inv_square_dist = (1 / d(i,j)^2)
                    DenseMatrix<DataType_> inv_square_dist(NodeDistanceInverse<>::value(_coordinates));

                    // Having inv_square_dist, mask it and invert each element to get square_dist
                    DenseMatrix<DataType_> square_dist(*inv_square_dist.copy());
                    ElementProduct<>::value(square_dist, _neighbours);
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
                    Scale<>::value(DataType_(1 / (_edge_length * _edge_length)), attractive_forces);

                    // Calculating inv_square_dist <- -inv_square_dist
                    Scale<DataType_>::value(DataType_(-1), inv_square_dist);

                    // Calculating diff = (diag(sum_vec_inv) - inv_square_dist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv->size(), *sum_vec_inv);
                    DenseMatrix<DataType_> diff(Sum<>::value(inv_square_dist, inv_diag));

                    // Calculating repulsive_forces = (_coordinates * diff)
                    DenseMatrix<DataType_> repulsive_forces(Product<>::value(_coordinates, diff));

                    // Calculating repulsive_forces = (k^2 * repulsive_forces)
                    Scale<>::value(DataType_(_edge_length * _edge_length), repulsive_forces);

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

        template <typename DataType_>
        class Implementation<DataType_, WeightedFruchtermanReingold>
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
                friend class Positions<DataType_, WeightedFruchtermanReingold>;

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
                    Scale<DataType_>::value(DataType_(-1), inv_square_dist);

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

        template <typename DataType_>
        class Implementation<DataType_, KamadaKawai>
        {
            private:
                /// position matrix - the coordinates of each node
                DenseMatrix<DataType_> & _coordinates;

                /// adjacence matrix - which nodes are neighbours?
                const DenseMatrix<bool> & _neighbours;

                /// a scalar factor which determinates how much influence the distance has while calculating the forces.
                const DataType_ _edge_length;

                /// parameter matrix for attractive forces
                DenseMatrix<DataType_> _attractive_force_parameter;

            public:
                friend class Positions<DataType_, KamadaKawai>;

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseMatrix<bool> & neighbours, const DataType_ edge_length) :
                    _coordinates(coordinates),
                    _neighbours(neighbours),
                    _edge_length(edge_length),
                    _attractive_force_parameter(neighbours.columns(), neighbours.rows(), 0)

                {
                    // Use adjacency matrix to calculate the _attractive_force_parameter matrix
                    for (typename Matrix<bool>::ConstElementIterator e(_neighbours.begin_elements()),
                        e_end(_neighbours.end_elements()) ; e != e_end ; ++e)
                    {
                        *e == 0 ? _attractive_force_parameter[e.row()][e.column()] = 2 * _attractive_force_parameter.columns() :
                        _attractive_force_parameter[e.row()][e.column()] = 1 ;
                    }

                    // Using Dijkstra to calculate the graph distance matrix represented by _attractive_force_parameter
                    Dijkstra<DataType_>::value(_attractive_force_parameter);

                    // Mul(_attractive_force_parameter, _attractive_force_parameter) * _edge_length^2
                    ElementProduct<>::value(_attractive_force_parameter, _attractive_force_parameter);
                    Scale<DataType_>::value(DataType_(_edge_length * _edge_length), _attractive_force_parameter);

                    // Calculate 1 / _attractive_force_parameter
                    ElementInverse<>::value(_attractive_force_parameter);
                }

                void init()
                {
                    //TODO: Implement KamadaKawai initialisation
                }


                DataType_ value(DataType_ & eps)
                {
                    // Calculate square_dist = d(i,j)^2
                    DenseMatrix<DataType_> square_dist(NodeDistance<>::value(_coordinates));

                    // Mul(square_dist, cost_matrix) -1
                    ElementProduct<>::value(square_dist, _attractive_force_parameter);
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

        template <typename DataType_>
        class Implementation<DataType_, WeightedKamadaKawai>
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
                friend class Positions<DataType_, WeightedKamadaKawai>;
                
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
