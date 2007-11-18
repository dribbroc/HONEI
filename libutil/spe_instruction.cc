/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libutil/assertion.hh>
#include <libutil/log.hh>
#include <libutil/spe_instruction.hh>
#include <libutil/spe_kernel.hh>

namespace honei
{
    SPEInstruction::SPEInstruction(const OpCode opcode, const unsigned size, const Operand a,
            const Operand b, const Operand c, const Operand d, const Operand e,
            const Operand f, const Operand g, const Operand h, const Operand i,
            const Operand j, const Operand k, const Operand l, const Operand m,
            const Operand n, const Operand o) :
        _kernel(0),
        _index(0)
    {
        _instruction.opcode = opcode;
        _instruction.size = size;
        _instruction.a = a;
        _instruction.b = b;
        _instruction.c = c;
        _instruction.d = d;
        _instruction.e = e;
        _instruction.f = f;
        _instruction.g = g;
        _instruction.h = h;
        _instruction.i = i;
        _instruction.j = j;
        _instruction.k = k;
        _instruction.l = l;
        _instruction.m = m;
        _instruction.n = n;
        _instruction.o = o;
    }

    SPEInstruction::~SPEInstruction()
    {
        wait();
    }

    void
    SPEInstruction::enqueue_with(SPEKernel * kernel) const
    {
        std::string msg("Enqueing SPEInstruction, opcode = " + stringify(_instruction.opcode) +
                ", size = " + stringify(_instruction.size) + "\n");
        msg += "a = " + stringify(_instruction.a.ea) + ", b = " + stringify(_instruction.b.ea) + "\n";
        msg += "c = " + stringify(_instruction.c.ea) + ", d = " + stringify(_instruction.d.ea) + "\n";
        msg += "e = " + stringify(_instruction.e.ea) + ", f = " + stringify(_instruction.f.ea) + "\n";
        msg += "g = " + stringify(_instruction.g.ea) + ", h = " + stringify(_instruction.h.ea) + "\n";
        msg += "i = " + stringify(_instruction.i.ea) + ", j = " + stringify(_instruction.j.ea) + "\n";
        msg += "k = " + stringify(_instruction.k.ea) + ", l = " + stringify(_instruction.l.ea) + "\n";
        msg += "m = " + stringify(_instruction.k.ea) + ", n = " + stringify(_instruction.l.ea) + "\n";
        msg += "o = " + stringify(_instruction.m.ea) + "\n";
        LOGMESSAGE(ll_minimal, msg);

        _kernel = kernel;
        _index = kernel->enqueue(_instruction);
    }

    void
    SPEInstruction::wait()
    {
        if (_kernel)
        {
            _kernel->wait(_index);
            _kernel = 0;
        }
    }

    const Operand
    SPEInstruction::empty = { 0 };
}
