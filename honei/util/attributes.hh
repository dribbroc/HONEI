/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
#ifndef HONEI_GUARD_ATTRIBUTES_HH
#define HONEI_GUARD_ATTRIBUTES_HH 1

#if defined (__GNUC__)
#  define HONEI_ALIGNED(x) __attribute__((__aligned__(x)))
#  define HONEI_INLINE __attribute__((__always_inline__))
#  define HONEI_PACKED __attribute__((packed))
#  define HONEI_THREAD_LOCAL static __thread
#  define HONEI_UNUSED __attribute__((unused))
#elif defined (DOXYGEN)
#  define HONEI_ALIGNED(x)
#  define HONEI_INLINE
#  define HONEI_PACKED
#  define HONEI_THREAD_LOCAL static
#  define HONEI_UNUSED
#else
#  error "Your compiler is not supported yet!"
#endif

#endif
