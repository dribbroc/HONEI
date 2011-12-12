/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef CUDA_GUARD_TRANSFER_HH
#define CUDA_GUARD_TRANSFER_HH 1

extern "C"
{
    void * cuda_malloc_host(unsigned long bytes);

    void * cuda_malloc(unsigned long bytes);

    void cuda_upload(void * src, void * target, unsigned long bytes);

    void cuda_download(void * src, void * target, unsigned long bytes);

    void cuda_free(void *src);

    void cuda_copy(void * src, void * dest, unsigned long bytes);

    void cuda_convert_float_double(void * src, void * dest, unsigned long bytes);

    void cuda_convert_double_float(void * src, void * dest, unsigned long bytes);

    void cuda_fill_zero(void * dest, unsigned long bytes);
}
#endif
