/* vim: set sw=4 sts=4 et nofoldenable syntax on */

/*
 * Copyright (c) 2007 Thorsten Deinert <thosten.deinert@uni-dortmund.de>
 *
 * This file is part of the Graph C++ library. LibGraph is free software;
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

#ifndef LIBGRAPH_GUARD_ABSTRACT_GRAPH
#define LIBGRAPH_GUARD_ABSTRACT_GRAPH 1

#include <tr1/memory>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/graph/node.hh>
#include <map>
#include <fstream>


namespace honei 
{
    template <typename DataType_> 
    class AbstractGraph
    {
    protected:
        DenseMatrix<DataType_> * _coordinates;
        SparseMatrix< DataType_ > * _edges;
        DenseVector<DataType_> * _node_weights;
    public:
        AbstractGraph() :
            _coordinates(0),
            _edges(0),
            _node_weights(0)
        {
        }

        virtual ~AbstractGraph()
        {
            delete _coordinates;
            delete _edges;
            delete _node_weights;
        }

        virtual inline DenseMatrix<DataType_> * coordinates()
        {
            return _coordinates;
        }

        virtual inline DenseVector<DataType_> * node_weights()
        {
            return _node_weights;
        }

        virtual inline SparseMatrix<DataType_> * edges()
        {
            return _edges;
        }

        virtual int node_count() = 0;
        
        virtual bool includes_timeslices()
        {
            return false;
        }
        
        virtual inline int slice_count()
        {
            return 0;
        }
                
        virtual inline int timeslice_index(int node_index)
        {
            return 0;
        }

        virtual inline bool same_timeslice(int index1, int index2)
        {
            return true;
        }
        
        void write_gml(const char filename[], bool include_coordinates = false)
        {
            std::ofstream fs(filename, std::ios_base::out);
            fs << "graph [\n";
            fs << "   directed 0\n";
            fs << "   timeslices " << (includes_timeslices() ? "1" : "0") << "\n";
            fs << "   coordinates " << (include_coordinates ? "1" : "0") << "\n";
            for (int i(0); i < _node_weights->size(); ++i)
            {
                fs << "    node [\n";
                fs << "        id " << i << "\n";
                fs << "        weight " << (*_node_weights)[i] << "\n";
                if (includes_timeslices())
                    fs << "        timeslice " << timeslice_index(i) << "\n";
                if (include_coordinates)
                {
                    fs << "        x " << (*_coordinates)(i, 0) << "\n";
                    fs << "        y " << (*_coordinates)(i, 1) << "\n";
                }
                fs << "    ]\n";
            }
            for (typename SparseMatrix<DataType_>::ConstElementIterator i(_edges->begin_elements()), end_i(_edges->end_elements());
                i != end_i; ++i)
                if (i.column() > i.row() && *i > 0)
                {
                        fs << "    edge [\n";
                        fs << "        source " << i.column() << "\n";
                        fs << "        target " << i.row() << "\n";
                        fs << "        weight " << *i << "\n";
                        fs << "    ]\n";
                }
            fs << "]\n";                
        fs.close();
        }
    };
}

#endif

