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

#include <honei/util/tags.hh>
#include <honei/util/exception.hh>
#include <honei/util/stringify.hh>

using namespace honei;

const std::string tags::CPU::name = "cpu";
const std::string tags::CPU::SSE::name = "sse";
const std::string tags::CPU::MultiCore::name = "mc";
const std::string tags::CPU::MultiCore::SSE::name = "mc-sse";
const std::string tags::Cell::name = "cell";
const std::string tags::GPU::name = "gpu";
const std::string tags::GPU::CUDA::name = "cuda";
const std::string tags::NONE::name = "none";

std::ostream & honei::operator<< (std::ostream & left, tags::TagValue value)
{
    do
    {
        switch (value)
        {
            case tags::tv_cpu:
                left << "CPU";
                continue;

            case tags::tv_cpu_multi_core:
                left << "CPU-Multicore";
                continue;

            case tags::tv_cell:
                left << "Cell";
                continue;

            case tags::tv_gpu:
                left << "GPU";
                continue;

            case tags::tv_gpu_cuda:
                left << "GPU CUDA";
                continue;

            case tags::tv_none:
                left << "NONE";
                continue;

            case tags::tv_fake:
                left << "Fake -- for test purpose only!";
                continue;
        }

        throw InternalError("Unexpected value for tags::TagValue '" + stringify(long(value)) + "'");
    } while (false);

    return left;
}

