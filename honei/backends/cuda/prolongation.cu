
/* vim: set sw=4 sts=4 et nofoldenable : */

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

#include <honei/backends/cuda/cuda_util.hh>

namespace honei
{
    namespace cuda
    {
        __global__ void prolongation_gpu(float * fine, float * coarse, int Nfine, int Ncoarse, int Mfine, int Mcoarse,
                int n1, int n2, int n3, int n4, int e1, int e2, int e3, int e4)
        {
            // index in the fine array (output index)
            int ifine = blockDim.x*blockIdx.x+threadIdx.x;
            // row index in the fine array  // TODO integer division
            int rfine = (int)floorf(ifine / (float)Mfine);
            // column index in the fine array
            int cfine = ifine - rfine*Mfine;
            // row index in the coarse array
            int rcoarse = (int)floorf(0.5f*rfine);
            // column index in the coarse array
            int ccoarse = (int)floorf(0.5f*cfine);
            // finally, index in the coarse array to start reading from
            // (bottom left node of the coarse cell in case of edge or inner points)
            int icoarse = rcoarse*Mcoarse + ccoarse;

            // compute odd/even information (==0 = even)
            int rodd = rfine & 1;
            int codd = cfine & 1;

            if (rodd==0)
            {
                if (codd==0)
                    // case 1: node is present in coarse and fine array
                    fine[ifine] = coarse[icoarse];
                else
                    // case 2: node on coarse edge, horizontal
                    fine[ifine] = 0.5f*(coarse[icoarse] + coarse[icoarse+1]);
            }
            else
            {
                if (codd==0)
                    // case 3: node on coarse edge, vertical
                    fine[ifine] = 0.5f*(coarse[icoarse] + coarse[icoarse+Mcoarse]);
                else
                    // case 4: inner node
                    fine[ifine] = 0.25f*(coarse[icoarse] + coarse[icoarse+1] + coarse[icoarse+Mcoarse] + coarse[icoarse+Mcoarse+1]);
            }
            // apply homogeneous Dirichlet BCs if necessary
            // quick coarse branch to minimise warp divergence
            if (rfine==0 || rfine==Mfine-1 || cfine==0 || cfine==Mfine-1)
            {
                if (ifine==0 && n1==2)
                    // bottom left node
                    fine[ifine] = 0.0f;
                else if (ifine==Mfine-1 && n2==2)
                    // bottom right node
                    fine[ifine] = 0.0f;
                else if (ifine==Nfine-1 && n3==2)
                    // top right node
                    fine[ifine] = 0.0f;
                else if (ifine==Nfine-Mfine && n4==2)
                    // top left node
                    fine[ifine] = 0.0f;
                else if (rfine==0 && cfine!=0 && cfine!=Mfine-1 && e1==2)
                    // bottom edge
                    fine[ifine] = 0.0f;
                else if (cfine==Mfine-1 && rfine!=0 && rfine!=Mfine-1 && e2==2)
                    // right edge
                    fine[ifine] = 0.0f;
                else if (rfine==Mfine-1 && cfine!=0 && cfine!=Mfine-1 && e3==2)
                    // top edge
                    fine[ifine] = 0.0f;
                else if (cfine==0 && rfine!=0 && rfine!=Mfine-1 && e4==2)
                    // left edge
                    fine[ifine] = 0.0f;
            }

        }
    }
}

extern "C" void cuda_prolongation_float(void * fine, unsigned long size_fine, void * coarse, unsigned long size_coarse,
        unsigned long * macroBorderMask, unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    unsigned long Nc = size_coarse;
    unsigned long Nf = size_fine;
    unsigned long Mc = (unsigned long)sqrt((double)Nc);
    unsigned long Mf = (unsigned long)sqrt((double)Nf);
    block.x  = blocksize;

    float * fine_gpu((float *)fine);
    float * coarse_gpu((float *)coarse);

    // directly apply BCs
    grid.x = (unsigned long)ceil(Nf/(double)block.x);
    honei::cuda::prolongation_gpu<<<grid,block>>>(fine_gpu, coarse_gpu, Nf, Nc, Mf, Mc,
            macroBorderMask[0], macroBorderMask[1], macroBorderMask[2], macroBorderMask[3],
            macroBorderMask[4], macroBorderMask[5], macroBorderMask[6], macroBorderMask[7]);
    CUDA_ERROR();
}
