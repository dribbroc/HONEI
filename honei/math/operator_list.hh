/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef MATH_GUARD_OPERATOR_LIST_HH
#define MATH_GUARD_OPERATOR_LIST_HH 1

#include <honei/util/private_implementation_pattern.hh>
#include <honei/math/operator.hh>
#include <vector>

namespace honei
{
    class OperatorList;
    template <> struct Implementation<OperatorList>
    {
        std::vector<Operator*> _ops;

        /// Constructor
        Implementation()
        {
        }

        /// Destructor
        ~Implementation()
        {
            for (unsigned long i(0) ; i < _ops.size() ; ++i)
                delete _ops.at(i);
        }
    };

    class OperatorList :
        public PrivateImplementationPattern<OperatorList, Shared>
    {

        public:
            OperatorList() :
                PrivateImplementationPattern<OperatorList, Shared>(new Implementation<OperatorList>())
            {
                CONTEXT("When constructing OperatorList:");
            }

            ~OperatorList()
            {
                CONTEXT("When destructing OperatorList:");
            }

            void push_back(Operator * op)
            {
                _imp->_ops.push_back(op);
            }

            void erase(unsigned long i)
            {
                _imp->_ops.erase(_imp->_ops.begin() + i);
            }

            void value()
            {
                CONTEXT("When running OperatorList::value():");
                for (unsigned long i(0) ; i < _imp->_ops.size() ; ++i)
                {
                    //std::cout << "CYCLE-SATUS: " << _imp->_ops.at(i)->to_string() << std::endl;
                    _imp->_ops.at(i)->value();
                }
            }

            Operator* operator[](unsigned long i)
            {
                return _imp->_ops.at(i);
            }

            unsigned long size()
            {
                return _imp->_ops.size();
            }
    };
}
#endif
