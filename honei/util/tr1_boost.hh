/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef HONEI_GUARD_TR1_BOOST_HH
#define HONEI_GUARD_TR1_BOOST_HH 1

#ifdef HONEI_BOOST
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/mem_fn.hpp>
#include <boost/type_traits.hpp>
#define HONEI_PLACEHOLDERS_1 _1
#define HONEI_PLACEHOLDERS_2 _2
#define HONEI_PLACEHOLDERS_3 _3
#define HONEI_PLACEHOLDERS_4 _4
#define HONEI_PLACEHOLDERS_5 _5
#define HONEI_PLACEHOLDERS_6 _6
using namespace boost;
#else
#include <tr1/memory>
#include <tr1/functional>
#include <tr1/type_traits>
#define HONEI_PLACEHOLDERS_1 std::tr1::placeholders::_1
#define HONEI_PLACEHOLDERS_2 std::tr1::placeholders::_2
#define HONEI_PLACEHOLDERS_3 std::tr1::placeholders::_3
#define HONEI_PLACEHOLDERS_4 std::tr1::placeholders::_4
#define HONEI_PLACEHOLDERS_5 std::tr1::placeholders::_5
#define HONEI_PLACEHOLDERS_6 std::tr1::placeholders::_6
using namespace std::tr1;
#endif

#endif
