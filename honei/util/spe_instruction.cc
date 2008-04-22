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

#include <honei/util/assertion.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/spe_instruction.hh>
#include <honei/util/spe_kernel.hh>

#include <tr1/functional>

namespace honei
{
    template <> struct Implementation<SPEInstruction>
    {
        typedef cell::Instruction Instruction;
        typedef cell::OpCode OpCode;
        typedef cell::Operand Operand;

        /// Our mutex.
        Mutex * const mutex;

        /// Our condition variable.
        ConditionVariable * const condition;

        /// Our are-we-finished variable.
        bool finished;

        /// Our instruction structure.
        Instruction instruction;

        /// Constructor.
        Implementation(const OpCode opcode, const unsigned size, const Operand & a,
                const Operand & b, const Operand & c, const Operand & d, const Operand & e,
                const Operand & f, const Operand & g, const Operand & h, const Operand & i,
                const Operand & j, const Operand & k, const Operand & l, const Operand & m,
                const Operand & n, const Operand & o) :
            mutex(new Mutex),
            condition(new ConditionVariable),
            finished(false)
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

        ~Implementation()
        {
            delete mutex;
            delete condition;
        }

        inline void wait()
        {
            Lock l(*mutex);

            while (! finished)
            {
                condition->wait(*mutex);
            }
        }
    };
}

using namespace honei;

SPEInstruction::SPEInstruction(const OpCode opcode, const unsigned size, const Operand & a,
        const Operand & b, const Operand & c, const Operand & d, const Operand & e,
        const Operand & f, const Operand & g, const Operand & h, const Operand & i,
        const Operand & j, const Operand & k, const Operand & l, const Operand & m,
        const Operand & n, const Operand & o) :
    PrivateImplementationPattern<SPEInstruction, Shared>(new Implementation<SPEInstruction>(
                opcode, size, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o))
{
}

SPEInstruction::~SPEInstruction()
{
}

SPEInstruction::Instruction &
SPEInstruction::instruction() const
{
    return _imp->instruction;
}

void
SPEInstruction::set_finished()
{
    Lock l(*_imp->mutex);

    _imp->finished = true;
    _imp->condition->signal();
}

void
SPEInstruction::wait() const
{
    _imp->wait();
}

bool
SPEInstruction::finished() const
{
    Lock l(*_imp->mutex);

    return _imp->finished;
}

const SPEInstruction::Operand
SPEInstruction::empty = { static_cast<void *>(0) };

template <typename DataType_>
SPEFrameworkInstruction<1, DataType_, cell::rtm_dma>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * elements, const unsigned size, const DataType_ scalar, const unsigned quantisation) :
    SPEInstruction(opcode, ((16384 / quantisation) / 16) * 16 * quantisation, elements),
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
        instruction.b.u = (size - skip) * sizeof(DataType_) / instruction.size;
        instruction.c.u = (size - skip) * sizeof(DataType_) % instruction.size;
        instruction.c.u &= ~0xF;
        instruction.c.u = ((instruction.c.u / quantisation) / 16) * quantisation * 16;

        _begin_transfers = skip;
        _end_transfers = ((instruction.b.u * instruction.size + instruction.c.u) / sizeof(DataType_)) + skip;

        if (sizeof(DataType_) == 4) // float
        {
            instruction.d.f = scalar;
        }
        else // double
        {
            instruction.d.d = scalar;
        }
    }

    if (0 == instruction.c.u)
    {
        if (instruction.b.u > 0)
        {
            instruction.c.u = instruction.size;
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
        const unsigned size, const DataType_ scalar) :
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
        instruction.c.u = (size - skip) * sizeof(DataType_) / instruction.size;
        instruction.d.u = (size - skip) * sizeof(DataType_) % instruction.size;
        instruction.d.u &= ~0xF;

        _begin_transfers = skip;
        _end_transfers = (instruction.c.u * instruction.size + instruction.d.u) / sizeof(DataType_) + skip;

        if (sizeof(DataType_) == 4) // float
        {
            instruction.e.f = scalar;
        }
        else // double
        {
            instruction.e.d = scalar;
        }
    }

    if (0 == instruction.d.u)
    {
        if (instruction.c.u > 0)
        {
            instruction.d.u = instruction.size;
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

template class SPEFrameworkInstruction<1, double, cell::rtm_dma>;

template class SPEFrameworkInstruction<1, double, cell::rtm_mail>;

template <typename DataType_>
SPEFrameworkInstruction<2, DataType_, cell::rtm_dma>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * a_elements,
        DataType_ * b_elements, const unsigned size, const DataType_ scalar) :
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
        unsigned skip((16 / sizeof(DataType_) - a_offset) % (16 / sizeof(DataType_)));
        instruction.b.u += (sizeof(DataType_) * skip); // Adjust SPU start for b respecting the elements calculated on PPU.
        unsigned b_offset((instruction.b.u & 0xF) / sizeof(DataType_)); // Alignment offset of b -> shuffle-factor on SPU.
        instruction.e.u = b_offset;

        // Align the address for SPU.
        instruction.a.u += (sizeof(DataType_) * skip);
        instruction.b.u += (sizeof(DataType_) * ((16 / sizeof(DataType_)) - b_offset));

        // Transmit the first __vector__ following the elements calculated on PPU (b_carry on SPU)
        Operand * dma_start(reinterpret_cast<Operand *>(instruction.b.ea));
        instruction.f = dma_start[-2];
        instruction.g = dma_start[-1];

        // Transmit scalar value.
        instruction.h = Operand(scalar);

        // Subtract PPU-calculated parts from size.
        instruction.c.u = (size - skip) / (16384 / sizeof(DataType_));
        instruction.d.u = (size - skip) % (16384 / sizeof(DataType_));
        instruction.d.u &= ~0xF;

        // Rest index dependent on offset and SPU part.
        _begin_transfers = skip;
        _end_transfers = (instruction.c.u * (16384 / sizeof(DataType_))) + instruction.d.u + skip;

        instruction.d.u *= sizeof(DataType_);
    }

    if (0 == instruction.d.u)
    {
        if (instruction.c.u > 0)
        {
            instruction.d.u = instruction.size;
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
template class SPEFrameworkInstruction<2, double, cell::rtm_dma>;

template <typename DataType_>
SPEFrameworkInstruction<2, DataType_, cell::rtm_mail>::SPEFrameworkInstruction(const OpCode opcode, DataType_ * result, DataType_ * a_elements,
        DataType_ * b_elements, const unsigned size, const DataType_ scalar) :
    SPEInstruction(opcode, 16384, result, a_elements, b_elements),
    _use_spe(true)
{
    Instruction & instruction(_imp->instruction);

    if (size < 5)
    {
        instruction.d.u = 0;
        instruction.e.u = 0;
        _begin_transfers = 0;
        _end_transfers = 0;
        _use_spe = false;
    }
    else
    {
        unsigned a_offset((instruction.b.u & 0xF) / sizeof(DataType_)); // Alignment offset of a -> elements calculated on PPU.
        unsigned skip((16 / sizeof(DataType_) - a_offset) % (16 / sizeof(DataType_)));
        instruction.c.u += (sizeof(DataType_) * skip); // Adjust SPU start for b respecting the elements calculated on PPU.
        unsigned b_offset((instruction.c.u & 0xF) / sizeof(DataType_)); // Alignment offset of b -> shuffle-factor on SPU.
        instruction.f.u = b_offset;

        // Align the address for SPU.
        instruction.b.u += (sizeof(DataType_) * skip);
        instruction.c.u += (sizeof(DataType_) * ((16 / sizeof(DataType_)) - b_offset));

        // Transmit the first __vector__ following the elements calculated on PPU (b_carry on SPU)
        Operand * dma_start(reinterpret_cast<Operand *>(instruction.c.ea));
        instruction.g = dma_start[-2];
        instruction.h = dma_start[-1];
        instruction.i = Operand(scalar);

        // Subtract PPU-calculated parts from size.
        instruction.d.u = (size - skip) / (16384 / sizeof(DataType_));
        instruction.e.u = (size - skip) % (16384 / sizeof(DataType_));
        instruction.e.u &= ~0xF;

        // Rest index dependent on offset and SPU part.
        _begin_transfers = skip;
        _end_transfers = (instruction.d.u * (16384 / sizeof(DataType_))) + instruction.e.u + skip;

        instruction.e.u *= sizeof(DataType_);
    }

    if (0 == instruction.e.u)
    {
        if (instruction.d.u > 0)
        {
            instruction.e.u = 16 * 1024;
        }
        else
        {
            _use_spe = false;
        }
    }
    else
    {
        ++instruction.d.u;
    }

}
template class SPEFrameworkInstruction<2, float, cell::rtm_mail>;
template class SPEFrameworkInstruction<2, double, cell::rtm_mail>;

template <typename DataType_>
SPEFrameworkInstruction<3, DataType_, cell::rtm_dma>::SPEFrameworkInstruction(const OpCode opcode,
        DataType_ * a_elements, DataType_ * b_elements, DataType_ * c_elements,
        const unsigned size, const DataType_ scalar) :
    SPEInstruction(opcode, 16384, a_elements, b_elements, c_elements), _use_spe(true)
{
    Instruction & instruction(_imp->instruction);

    if (size < 5)
    {
        instruction.d.u = 0;
        instruction.e.u = 0;
        _begin_transfers = 0;
        _end_transfers = 0;
        _use_spe = false;
    }
    else
    {
        unsigned a_offset((instruction.a.u & 0xF) / sizeof(DataType_)); // Alignment offset of a -> elements calculated on PPU.
        unsigned skip((16 / sizeof(DataType_) - a_offset) % (16 / sizeof(DataType_)));
        instruction.b.u += (sizeof(DataType_) * skip); // Adjust SPU start for b respecting the elements calculated on PPU.
        instruction.c.u += (sizeof(DataType_) * skip); // Adjust SPU start for c respecting the elements calculated on PPU.
        unsigned b_offset((instruction.b.u & 0xF) / sizeof(DataType_)); // Alignment offset of b -> shuffle-factor on SPU.
        unsigned c_offset((instruction.c.u & 0xF) / sizeof(DataType_)); // Alignment offset of c -> shuffle-factor on SPU.
        instruction.f.u = b_offset;
        instruction.g.u = c_offset;

        // Align the address for SPU.
        instruction.a.u += (sizeof(DataType_) * skip);
        instruction.b.u += (sizeof(DataType_) * ((16 / sizeof(DataType_)) - b_offset));
        instruction.c.u += (sizeof(DataType_) * ((16 / sizeof(DataType_)) - c_offset));

        // Transmit the first __vector__ following the elements calculated on PPU (b_carry and c_carry on SPU)
        DataType_ * b_dma_start = reinterpret_cast<DataType_ *>(instruction.b.u);
        DataType_ * c_dma_start = reinterpret_cast<DataType_ *>(instruction.c.u);
        if (sizeof(DataType_) == 4) // float
        {
            instruction.h.f = *(b_dma_start - 4);
            instruction.i.f = *(b_dma_start - 3);
            instruction.j.f = *(b_dma_start - 2);
            instruction.k.f = *(b_dma_start - 1);
            instruction.l.f = *(c_dma_start - 4);
            instruction.m.f = *(c_dma_start - 3);
            instruction.n.f = *(c_dma_start - 2);
            instruction.o.f = *(c_dma_start - 1);
            instruction.p.f = scalar;
        }
        else // double
        {
            instruction.h.d = *(b_dma_start - 2);
            instruction.i.d = *(b_dma_start - 1);
            instruction.l.d = *(c_dma_start - 2);
            instruction.m.d = *(c_dma_start - 1);
            instruction.p.d = scalar;
        }

        // Subtract PPU-calculated parts from size.
        instruction.d.u = (size - skip) / (16384 / sizeof(DataType_));
        instruction.e.u = (size - skip) % (16384 / sizeof(DataType_));
        instruction.e.u &= ~0xF;

        // Rest index dependent on offset and SPU part.
        _begin_transfers = skip;
        _end_transfers = (instruction.d.u * (16384 / sizeof(DataType_))) + instruction.e.u + skip;

        instruction.e.u *= sizeof(DataType_);
    }

    if (0 == instruction.e.u)
    {
        if (instruction.d.u > 0)
        {
            instruction.e.u = instruction.size;
        }
        else
        {
            _use_spe = false;
        }
    }
    else
    {
        ++instruction.d.u;
    }

}
template class SPEFrameworkInstruction<3, float, cell::rtm_dma>;
template class SPEFrameworkInstruction<3, double, cell::rtm_dma>;

namespace honei
{
    template <> struct Implementation<SPEInstructionQueue>
    {
        /// Our list of instructions.
        std::list<SPEInstruction> instructions;
    };
}

SPEInstructionQueue::SPEInstructionQueue() :
    PrivateImplementationPattern<SPEInstructionQueue, Single>(new Implementation<SPEInstructionQueue>)
{
}

SPEInstructionQueue::~SPEInstructionQueue()
{
    wait();
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

unsigned long
SPEInstructionQueue::size()
{
    return _imp->instructions.size();
}

bool
SPEInstructionQueue::finished() const
{
    Lock l(*_imp->instructions.back()._imp->mutex);

    return _imp->instructions.back()._imp->finished;
}
