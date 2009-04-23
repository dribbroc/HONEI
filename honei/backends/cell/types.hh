/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/util/attributes.hh>

namespace honei
{
    namespace cell
    {

#if defined(__PPU__)

        struct LocalStoreAddress
        {
            unsigned value;
        };

        typedef void * EffectiveAddress;

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

        template <> union Pointer<unsigned>
        {
            void * untyped;
            unsigned adjustable;
            unsigned * typed;
            vector unsigned * vectorised;
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

        template <> union Subscriptable<double>
        {
            vector double value;
            double array[2];
        };


        template <> union Subscriptable<unsigned long>
        {
            vector unsigned long value;
            unsigned long array[4];
        };

        template <> union Subscriptable<unsigned>
        {
            vector unsigned value;
            unsigned array[4];
        };

        template <typename T_> union Vector;

        template <> union Vector<float>
        {
            vector float vf;
            vector unsigned int vui;
        };

        template <> union Vector<double>
        {
            vector double vd;
            vector unsigned int vui;
        };

        template <> union Vector<unsigned>
        {
            vector unsigned int vui;
            vector float vf;
        };

        template <typename T_> union MailableResult;

        template <> union MailableResult<float>
        {
            float value;
            unsigned int mail;
        };

        template <> union MailableResult<double>
        {
            double value;
            unsigned int mail;
        };
#else
# error "You should not include this header in anything but Cell source code!"
#endif

        struct HONEI_PACKED HONEI_ALIGNED(8) ListElement
        {
            unsigned stall_and_notify : 1;
            unsigned reserved : 15;
            unsigned transfer_size : 16;
            unsigned effective_address_low;
        };

        template <typename T_> union Pointer;

        template <> union Pointer<ListElement>
        {
            void * untyped;
            ListElement le;
        };

    }


}
#endif
