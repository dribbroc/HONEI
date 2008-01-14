/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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

SPEInstruction::Instruction &
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

template <typename DataType_>
SPEFrameworkInstruction<1, DataType_, cell::rtm_dma>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * elements, const unsigned size) :
    SPEInstruction(opcode, 16384, elements),
    _use_spe(true)
{
    Instruction & instruction(_imp->instruction);
    unsigned offset(instruction.a.u & 0xF), skip(0);

    if (offset > 0)
    {
        instruction.a.u += 16 - offset; // Align the address.
        skip = (16 - offset) / sizeof(DataType_);
    }

    if (size < 5)
    {
        offset = 0;
        instruction.b.u = 0;
        instruction.c.u = 0;
        _begin_transfers = 0;
        _end_transfers = 0;
    }
    else
    {
        instruction.b.u = (size - skip) / (16384 / sizeof(DataType_));
        instruction.c.u = (size - skip) % (16384 / sizeof(DataType_));
        instruction.c.u &= ~0xF;

        _begin_transfers = skip;
        _end_transfers = instruction.b.u * 4096 + instruction.c.u + skip;
        instruction.c.u *= 4;
    }

    if (0 == instruction.c.u)
    {
        if (instruction.b.u > 0)
        {
            instruction.c.u = 16 * 1024;
        }
        else
        {
            _use_spe = false;
        }
    }
    else
    {
        ++instruction.b.u;
    }
}

template <typename DataType_>
SPEFrameworkInstruction<1, DataType_, cell::rtm_mail>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * result, DataType_ * elements,
        const unsigned size) :
    SPEInstruction(opcode, 16384, result, elements),
    _use_spe(true)
{
    Instruction & instruction(_imp->instruction);
    unsigned offset(instruction.b.u & 0xF), skip(0);

    if (offset > 0)
    {
        instruction.b.u += 16 - offset; // Align the address.
        skip = (16 - offset) / sizeof(DataType_);
    }

    if (size < 5)
    {
        offset = 0;
        instruction.c.u = 0;
        instruction.d.u = 0;
        _begin_transfers = 0;
        _end_transfers = 0;
    }
    else
    {
        instruction.c.u = (size - skip) / (16384 / sizeof(DataType_));
        instruction.d.u = (size - skip) % (16384 / sizeof(DataType_));
        instruction.d.u &= ~0xF;

        _begin_transfers = skip;
        _end_transfers = instruction.c.u * 4096 + instruction.d.u + skip;
        instruction.d.u *= 4;
    }

    if (0 == instruction.d.u)
    {
        if (instruction.c.u > 0)
        {
            instruction.d.u = 16 * 1024;
        }
        else
        {
            _use_spe = false;
        }
    }
    else
    {
        ++instruction.c.u;
    }
}
template class SPEFrameworkInstruction<1, float, cell::rtm_dma>;

template class SPEFrameworkInstruction<1, float, cell::rtm_mail>;

template <typename DataType_>
SPEFrameworkInstruction<2, DataType_, cell::rtm_dma>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * a_elements,
        DataType_ * b_elements, const unsigned size) :
    SPEInstruction(opcode, 16384, a_elements, b_elements),
    _use_spe(true)
{
    Instruction & instruction(_imp->instruction);

    if (size < 5)
    {
        instruction.c.u = 0;
        instruction.d.u = 0;
        _begin_transfers = 0;
        _end_transfers = 0;
        _use_spe = false;
    }
    else
    {
        unsigned a_offset((instruction.a.u & 0xF) / sizeof(DataType_)); // Alignment offset of a -> elements calculated on PPU.
        unsigned a_inv_offset((4 - a_offset) % 4);
        instruction.b.u += (4 * a_inv_offset); // Adjust SPU start for b respecting the elements calculated on PPU.
        unsigned b_offset((instruction.b.u & 0xF) / sizeof(DataType_)); // Alignment offset of b -> shuffle-factor on SPU.
        instruction.e.u = b_offset;

        // Transmit the first __vector__ following the elements calculated on PPU (b_carry on SPU)
        DataType_ carries[4] = { b_elements[a_inv_offset], b_elements[a_inv_offset + 1],
            b_elements[a_inv_offset + 2], b_elements[a_inv_offset + 3] };

        for (unsigned i(0) ; i < b_offset ; i++)
        {
            carries[3] = carries[2];
            carries[2] = carries[1];
            carries[1] = carries[0];
            carries[0] = DataType_(0);
        }

        instruction.f.f = carries[0];
        instruction.g.f = carries[1];
        instruction.h.f = carries[2];
        instruction.i.f = carries[3];

        // Align the address for SPU.
        instruction.a.u += (4 * a_inv_offset);
        instruction.b.u += (4 * (4 - b_offset));

        //Subtract PPU-calculated parts from size.
        instruction.c.u = (size - a_inv_offset) / (1024 * 4);
        instruction.d.u = (size - a_inv_offset) % (1024 * 4);
        instruction.d.u &= ~0xF;

        // Rest index dependent on offset and SPU part.
        _begin_transfers = a_inv_offset;
        _end_transfers = (instruction.c.u * 4096) + instruction.d.u + a_inv_offset;

        instruction.d.u *= 4;
    }

    if (0 == instruction.d.u)
    {
        if (instruction.c.u > 0)
        {
            instruction.d.u = 16 * 1024;
        }
        else
        {
            _use_spe = false;
        }
    }
    else
    {
        ++instruction.c.u;
    }

}
template class SPEFrameworkInstruction<2, float, cell::rtm_dma>;

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
