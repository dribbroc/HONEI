/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
#include <libla/matrix_difference.hh>
#include <libla/matrix_error.hh>
#include <libla/matrix_element_inverse.hh>
#include <libla/matrix_elementwise_product.hh>
#include <libla/matrix_product.hh>
#include <libla/matrix_row_sum_vector.hh>
#include <libla/matrix_sum.hh>
#include <libla/scalar_matrix_product.hh>
#include <libla/scalar_matrix_sum.hh>
#include <libla/vector.hh>

#include <iostream>
#include <tr1/memory>
#include <cmath>

/**
 * \file
 *
 * Interface declarations for position calculation classes.
 *
 * \ingroup grplibgraph
 **/
 namespace pg512
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
            }

            void update(DataType_ eps, int timeout)
            {
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
                    /// \todo Move these checks to Position's constructor?
                    if (coordinates.columns() != neighbours.columns())
                        throw MatrixColumnsDoNotMatch(neighbours.columns(), coordinates.columns());

                    if (! neighbours.square())
                        throw MatrixIsNotSquare(neighbours.rows(), neighbours.columns());

#if 0
                    \todo Implement Trace
                    if (Trace<>::value(neighbours) > false)
                        throw MatrixHasNonVanishingTrace;

                    if (edge_length <= std::numeric_limits<DataType_>::epsilon())
                        \todo Implement GraphError
                        throw GraphError("Edge length must be positive");
#endif
                }

                DataType_ value(const DataType_ &)
                {
                    // Calculate InvSquaDist = (1 / d(i,j))
                    DenseMatrix<DataType_> inv_square_dist(NodeDistanceInverse<>::value(_coordinates));

                    // Having InvSquaDist, mask it and invert each element to get SquaDist
                    DenseMatrix<DataType_> square_dist(*inv_square_dist.copy());
                    MatrixElementwiseProduct<>::value(square_dist, _neighbours);
                    MatrixElementInverse<>::value(square_dist);

                    // Calculate the single diagonal matrix containing the row sum vector of SquaDist
                    DenseVector<DataType_> & sum_vec(*MatrixRowSumVector<DataType_>::value(square_dist));

                    // Calculate the same stuff for InvSquaDist
                    DenseVector<DataType_> & sum_vec_inv(*MatrixRowSumVector<DataType_>::value(inv_square_dist));

                    // Calculating SquaDist = (SquaDist - diag(SumVec))
                    BandedMatrix<DataType_> diag(sum_vec.size(), &sum_vec);
                    MatrixDifference<>::value(diag, square_dist); /// \todo Switch to MD::value(Dense, const Banded)
                    ScalarMatrixProduct<>::value(DataType_(-1), square_dist);

                    // Calculating FAttr = (p * SquaDist)
                    DenseMatrix<DataType_> attractive_forces(MatrixProduct<>::value(_coordinates, square_dist));

                    // Calculating FAttr = (1/k^2 * FAttr)
                    ScalarMatrixProduct<>::value(DataType_(1 / (_edge_length * _edge_length)), attractive_forces);

                    // Calculating inv_square_dist <- -inv_square_dist
                    ScalarMatrixProduct<DataType_>::value(DataType_(-1), inv_square_dist);

                    // Calculating Distance_Neg = (diag(SumVecInv) - InvSquaDist)
                    BandedMatrix<DataType_> inv_diag(sum_vec_inv.size(), &sum_vec_inv);
                    MatrixSum<>::value(inv_square_dist, diag);

                    // Calculating FRep = (p * Distance_neg)
                    DenseMatrix<DataType_> repulsive_forces(MatrixProduct<>::value(_coordinates, inv_diag));

                    // Calculating FRep = (k^2 * FRep)
                    ScalarMatrixProduct<>::value(_edge_length * _edge_length, repulsive_forces);

                    // Calculate the maximal force and the result forces
                    MatrixSum<>::value(attractive_forces, repulsive_forces);
                    DataType_ result(0);
                    DenseVector<DataType_> resulting_forces(_coordinates.columns(), DataType_(0)); /// \todo check
                    for (typename DenseVector<DataType_>::ElementIterator i(resulting_forces.begin_elements()), i_end(resulting_forces.end_elements()) ;
                        i != i_end ; ++i)
                    {
                        *i = VectorNorm<DataType_, vnt_l_two, true>::value(attractive_forces.column(i.index()));
                        result = std::max(result, *i);
                    }

                    // Calculate the new positions by using result forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_coordinates.begin_elements()),
                        e_end(_coordinates.end_elements()), k(attractive_forces.begin_elements()) ; e != e_end ; ++e, ++k)
                    {
                        *e = *e + 1 / (10 * resulting_forces[e.column()]) * *k;
                    }

                    return result;
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
                DataType_ _edge_length;

            public:
                friend class Positions<DataType_, KamadaKawai>;

                Implementation(DenseMatrix<DataType_> & coordinates, const DenseMatrix<bool> & neighbours, DataType_ edge_length) :
                    _coordinates(coordinates),
                    _neighbours(neighbours),
                    _edge_length(edge_length)
                {
                    /// \todo See comment in FruchtermanReingold-Implementation's
                    if (coordinates.columns() != neighbours.columns())
                        throw MatrixColumnsDoNotMatch(neighbours.columns(), coordinates.columns());

                    if (neighbours.columns() != neighbours.rows())
                        throw MatrixColumnsDoNotMatch(neighbours.rows(), neighbours.columns());
                }

                DataType_ value(DataType_ & eps)
                {
                    /// \todo Apply coding conventions.
                    // Calculate SquaDist = d(i,j)
                    DenseMatrix<DataType_> SquaDist = NodeDistance<>::value(_coordinates);

                    // Use adjacency matrix to calculate the cost matrix
                    DenseMatrix<DataType_> CostMat(_neighbours.columns(), _neighbours.columns(), 0);
                    int ix=0;
                    int jx=0;
                    for (typename Matrix<bool>::ConstElementIterator e(_neighbours.begin_elements()),
                        e_end(_neighbours.end_elements()) ; e != e_end ; ++e)
                    {
                        *e == 0 ? CostMat[ix][jx] = 2 * CostMat.columns() : CostMat[ix][jx] = 1 ;
                        jx = ((jx +1) % CostMat.columns());
                        jx == 0 ? ix++ : ix = ix;
                    }


                    // Calculate the cost matrix by Dijkstra
                    DenseMatrix<DataType_> newCostMat = Dijkstra<DataType_>::value(*CostMat.copy());

                    // Calculate InvGraphDist = Mul(newCostMat,newCostMat) * k^2
                    DenseMatrix<DataType_> InvGraphDist = MatrixElementwiseProduct<>::value (newCostMat, newCostMat);
                    InvGraphDist = ScalarMatrixProduct<DataType_>::value((_edge_length * _edge_length),InvGraphDist);

                    // Calculate InvGraphDist = 1 / InvGraphDist
                    InvGraphDist = MatrixElementInverse<>::value(InvGraphDist);

                    // Calculate SquaDist_InvGD = Mul(InvGraphDist,SquaDist) -1
                    DenseMatrix<DataType_> SquaDist_InvGD = ScalarMatrixSum<>::value (-1, MatrixElementwiseProduct<>::value (InvGraphDist, SquaDist));

                    // Calculate the row sum vector of SquaDist_InvGD, SumVec
                    std::tr1::shared_ptr<DenseVector<DataType_> > SumVec(MatrixRowSumVector<>::value(SquaDist_InvGD));

                    // Calculate SquaDist_InvGD = SquaDist_InvGD - diag(SumVec)
                    for (int i = 0; i < SquaDist_InvGD.columns(); i++)
                    {
                        SquaDist_InvGD[i][i] -= (*SumVec)[i];
                    }

                    // Calculate Forces: fKK = p * (SquaDist_InvGD - diag(SumVec))
                    DenseMatrix<DataType_> fKK =  MatrixProduct<>::value(*_coordinates.copy(), SquaDist_InvGD);

                    // Calculate the result forces, the maximal force and the node with maximal force
                    DenseVector<DataType_> fRes(_coordinates.columns(), 0, 0, 1);
                    DataType_ maxForce = 0;
                    int maxNode = 0;
                    for (int i = 0; i < fKK.columns(); ++i)
                    {
                        DenseVector<DataType_> column = fKK.column(i);
                        DataType_ f = VectorNorm<DataType_, vnt_l_two, true>::value(column);
                        fRes[i] = f;
                        maxForce < f ? maxNode = i, maxForce = f : maxNode = maxNode;
                    }

                    jx = 0;
                    ix = 0;

                    // Calculate the new positions by using result forces
                    for (typename MutableMatrix<DataType_>::ElementIterator e(_coordinates.begin_elements()),
                        e_end(_coordinates.end_elements()) ; e != e_end ; ++e)
                    {
                        (jx == maxNode) && (maxForce > eps)? *e = *e + 1 / (10 * fRes[maxNode]) * fKK[ix][maxNode]: *e = *e ;
                        jx = ((jx +1) % _coordinates.columns());
                        jx == 0 ? ix++ : ix = ix;
                    }
                    return maxForce;

                }
        };

    }
}

#endif
