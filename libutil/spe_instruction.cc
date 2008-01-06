/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <tr1/functional>

using namespace honei;

struct SPEInstruction::Implementation
{
    /// Our kernel.
    SPEKernel * kernel;

    /// Our index.
    unsigned index;

    /// Our instruction structure.
    Instruction instruction;

    /// Constructor.
    Implementation(const OpCode opcode, const unsigned size, const Operand & a,
            const Operand & b, const Operand & c, const Operand & d, const Operand & e,
            const Operand & f, const Operand & g, const Operand & h, const Operand & i,
            const Operand & j, const Operand & k, const Operand & l, const Operand & m,
            const Operand & n, const Operand & o) :
        kernel(0),
        index(0)
    {

        instruction.opcode = opcode;
        instruction.size = size;
        instruction.a = a;
        instruction.b = b;
        instruction.c = c;
        instruction.d = d;
        instruction.e = e;
        instruction.f = f;
        instruction.g = g;
        instruction.h = h;
        instruction.i = i;
        instruction.j = j;
        instruction.k = k;
        instruction.l = l;
        instruction.m = m;
        instruction.n = n;
        instruction.o = o;
    }

    inline void wait()
    {
        if (kernel)
        {
            kernel->wait(index);

            delete kernel;
            kernel = 0;
        }
    }

    inline bool finished()
    {
        if (kernel)
            return kernel->finished(index);
        else
            return false;
    }
};

SPEInstruction::SPEInstruction(const OpCode opcode, const unsigned size, const Operand & a,
        const Operand & b, const Operand & c, const Operand & d, const Operand & e,
        const Operand & f, const Operand & g, const Operand & h, const Operand & i,
        const Operand & j, const Operand & k, const Operand & l, const Operand & m,
        const Operand & n, const Operand & o) :
    _imp(new Implementation(opcode, size, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o))
{
}

SPEInstruction::~SPEInstruction()
{
    _imp->wait();
}

void
SPEInstruction::_enqueue_with(SPEKernel * kernel) const
{
    LOGMESSAGE(ll_minimal, std::string("Enqueing SPEInstruction, opcode = " + stringify(_imp->instruction.opcode)
                + ", size = " + stringify(_imp->instruction.size) + "\na = " + stringify(_imp->instruction.a.ea)
                + ", b = " + stringify(_imp->instruction.b.ea) + "\nc = " + stringify(_imp->instruction.c.ea)
                + ", d = " + stringify(_imp->instruction.d.ea) + "\ne = " + stringify(_imp->instruction.e.ea)
                + ", f = " + stringify(_imp->instruction.f.ea) + "\ng = " + stringify(_imp->instruction.g.ea)
                + ", h = " + stringify(_imp->instruction.h.ea) + "\ni = " + stringify(_imp->instruction.i.ea)
                + ", j = " + stringify(_imp->instruction.j.ea) + "\nk = " + stringify(_imp->instruction.k.ea)
                + ", l = " + stringify(_imp->instruction.l.ea) + "\nm = " + stringify(_imp->instruction.k.ea)
                + ", n = " + stringify(_imp->instruction.l.ea) + "\no = " + stringify(_imp->instruction.m.ea)
                + "\n"));

    _imp->kernel = new SPEKernel(*kernel);
    _imp->index = kernel->enqueue(*this);
}

SPEInstruction::Instruction
SPEInstruction::instruction() const
{
    return _imp->instruction;
}

void
SPEInstruction::wait() const
{
    _imp->wait();
}

bool
SPEInstruction::finished() const
{
    return _imp->finished();
}

const SPEInstruction::Operand
SPEInstruction::empty = { static_cast<void *>(0) };

struct SPEInstructionQueue::Implementation
{
    /// Our list of instructions.
    std::list<SPEInstruction> instructions;
};

SPEInstructionQueue::SPEInstructionQueue() :
    _imp(new Implementation)
{
}

SPEInstructionQueue::~SPEInstructionQueue()
{
    wait();
    delete _imp;
}

void
SPEInstructionQueue::push_back(const SPEInstruction & instruction)
{
    if (_imp->instructions.size() > 0)
    {
        if (instruction.instruction().opcode != (*(_imp->instructions.begin())).instruction().opcode)
            throw InternalError("SPEInstructionQueue: Inserting different opcodes in one Queue.");
    }
    if (_imp->instructions.size() > 7)
            throw InternalError("SPEInstructionQueue: InstructionQueue size is limited to 8 elements due to restrictions in the spu programm.");

    _imp->instructions.push_back(instruction);
}

void
SPEInstructionQueue::wait() const
{
    std::for_each(_imp->instructions.begin(), _imp->instructions.end(),
            std::tr1::bind(std::tr1::mem_fn(&SPEInstruction::wait), std::tr1::placeholders::_1));
}

std::list<SPEInstruction>::iterator
SPEInstructionQueue::begin()
{
    return _imp->instructions.begin();
}

std::list<SPEInstruction>::iterator
SPEInstructionQueue::end()
{
    return _imp->instructions.end();
}

void
SPEInstructionQueue::_enqueue_with(SPEKernel * kernel)
{
    if (_imp->instructions.size() > 0)
    {
        for (std::list<SPEInstruction>::iterator i_it(_imp->instructions.begin()), i_end(_imp->instructions.end()) ;
                i_it != i_end ; i_it++)
        {
            LOGMESSAGE(ll_minimal, std::string("(SPEInstructionQueue) Enqueing SPEInstruction, opcode = " + stringify((*i_it)._imp->instruction.opcode)
                        + ", size = " + stringify((*i_it)._imp->instruction.size) + "\na = " + stringify((*i_it)._imp->instruction.a.ea)
                        + ", b = " + stringify((*i_it)._imp->instruction.b.ea) + "\nc = " + stringify((*i_it)._imp->instruction.c.ea)
                        + ", d = " + stringify((*i_it)._imp->instruction.d.ea) + "\ne = " + stringify((*i_it)._imp->instruction.e.ea)
                        + ", f = " + stringify((*i_it)._imp->instruction.f.ea) + "\ng = " + stringify((*i_it)._imp->instruction.g.ea)
                        + ", h = " + stringify((*i_it)._imp->instruction.h.ea) + "\ni = " + stringify((*i_it)._imp->instruction.i.ea)
                        + ", j = " + stringify((*i_it)._imp->instruction.j.ea) + "\nk = " + stringify((*i_it)._imp->instruction.k.ea)
                        + ", l = " + stringify((*i_it)._imp->instruction.l.ea) + "\nm = " + stringify((*i_it)._imp->instruction.k.ea)
                        + ", n = " + stringify((*i_it)._imp->instruction.l.ea) + "\no = " + stringify((*i_it)._imp->instruction.m.ea)
                        + "\n"));

            (*i_it)._imp->kernel = new SPEKernel(*kernel);
            (*i_it)._imp->index = kernel->enqueue_queue_element(*i_it);
        }
        (*_imp->instructions.begin())._imp->kernel->run_queue();
    }
}

unsigned long
SPEInstructionQueue::size()
{
    return _imp->instructions.size();
}

bool
SPEInstructionQueue::finished() const
{
    return _imp->instructions.back()._imp->finished();
}

struct SPEInstructionStream::Implementation
{
    /// Our list of instructions.
    std::vector<SPEInstruction> instructions;

    /// Are we already enqueued with a kernel ?
    bool enqueued;

    /// The kernel we stream our instructions to.
    SPEKernel * kernel;
};

SPEInstructionStream::SPEInstructionStream(const SPEInstruction & instruction) :
    _imp(new Implementation)
{
    _imp->instructions.push_back(instruction);
    _imp->kernel = 0;
    _imp->enqueued = false;
}

SPEInstructionStream::~SPEInstructionStream()
{
    wait();
    delete _imp->kernel;
    delete _imp;
}

void
SPEInstructionStream::input(const SPEInstruction & instruction)
{
    if (instruction.instruction().opcode != (*(_imp->instructions.begin())).instruction().opcode)
        throw InternalError("SPEInstructionStream: Inserting different opcodes in one Stream.");

    _imp->instructions.push_back(instruction);
    if (_imp->enqueued)
    {
        _imp->instructions.back()._enqueue_with(_imp->kernel);
    }
}

void
SPEInstructionStream::wait() const
{
    std::for_each(_imp->instructions.begin(), _imp->instructions.end(),
            std::tr1::bind(std::tr1::mem_fn(&SPEInstruction::wait), std::tr1::placeholders::_1));
}

std::vector<SPEInstruction>::iterator
SPEInstructionStream::begin()
{
    return _imp->instructions.begin();
}

std::vector<SPEInstruction>::iterator
SPEInstructionStream::end()
{
    return _imp->instructions.end();
}

void
SPEInstructionStream::_enqueue_with(SPEKernel * kernel)
{
    if (_imp->enqueued || _imp->kernel)
        throw InternalError("SPEInstructionStream: You cannot enqueue one stream with any kernels twice.");

    _imp->kernel = new SPEKernel(*kernel);
    _imp->enqueued = true;
    for (std::vector<SPEInstruction>::iterator i_it(_imp->instructions.begin()), i_end(_imp->instructions.end()) ;
            i_it != i_end ; i_it++)
    {
        LOGMESSAGE(ll_minimal, std::string("(SPEInstructionStream) Enqueing SPEInstruction, opcode = " + stringify((*i_it)._imp->instruction.opcode)
                    + ", size = " + stringify((*i_it)._imp->instruction.size) + "\na = " + stringify((*i_it)._imp->instruction.a.ea)
                    + ", b = " + stringify((*i_it)._imp->instruction.b.ea) + "\nc = " + stringify((*i_it)._imp->instruction.c.ea)
                    + ", d = " + stringify((*i_it)._imp->instruction.d.ea) + "\ne = " + stringify((*i_it)._imp->instruction.e.ea)
                    + ", f = " + stringify((*i_it)._imp->instruction.f.ea) + "\ng = " + stringify((*i_it)._imp->instruction.g.ea)
                    + ", h = " + stringify((*i_it)._imp->instruction.h.ea) + "\ni = " + stringify((*i_it)._imp->instruction.i.ea)
                    + ", j = " + stringify((*i_it)._imp->instruction.j.ea) + "\nk = " + stringify((*i_it)._imp->instruction.k.ea)
                    + ", l = " + stringify((*i_it)._imp->instruction.l.ea) + "\nm = " + stringify((*i_it)._imp->instruction.k.ea)
                    + ", n = " + stringify((*i_it)._imp->instruction.l.ea) + "\no = " + stringify((*i_it)._imp->instruction.m.ea)
                    + "\n"));

        (*i_it)._imp->kernel = new SPEKernel(*(_imp->kernel));
        (*i_it)._imp->index = kernel->enqueue(*i_it);
    }
}

unsigned long
SPEInstructionStream::size()
{
    return _imp->instructions.size();
}

bool
SPEInstructionStream::finished() const
{
    return _imp->instructions.back()._imp->finished();
}
