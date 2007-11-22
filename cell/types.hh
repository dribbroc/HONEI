/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef CELL_GUARD_TYPES_HH
#define CELL_GUARD_TYPES_HH 1

namespace honei
{
    namespace cell
    {

#if defined(__PPU__)

        typedef void * EffectiveAddress;

        typedef unsigned int LocalStoreAddress;

#elif defined(__SPU__)

        typedef unsigned long long EffectiveAddress;

        typedef void * LocalStoreAddress;

        template <typename T_> union Pointer;

        template <> union Pointer<float>
        {
            void * untyped;
            unsigned adjustable;
            float * typed;
            vector float * vectorised;
        };

        template <> union Pointer<double>
        {
            void * untyped;
            unsigned adjustable;
            double * typed;
            vector double * vectorised;
        };

        template <> union Pointer<unsigned long long>
        {
            void * untyped;
            unsigned adjustable;
            unsigned long long * typed;
            vector unsigned long long * vectorised;
        };

        template <typename T_> union Subscriptable;

        template <> union Subscriptable<float>
        {
            vector float value;
            float array[4];
        };

        template <> union Subscriptable<unsigned long>
        {
            vector unsigned long value;
            unsigned long array[4];
        };

        template <typename T_> union MailableResult;

        template <> union MailableResult<float>
        {
            float value;
            unsigned int mail;
        };

#else
# error "You should not include this header in anything but Cell source code!"
#endif

    }
}

#endif
