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
        __global__ void restriction_gpu(float * fine, float * coarse, int Nc, int Nf, int Mc, int Mf,
                int n1, int n2, int n3, int n4, int e1, int e2, int e3, int e4)
        {
            // index in the coarse array (output index)
            int ic = blockDim.x*blockIdx.x+threadIdx.x;
            // row in the coarse array
            int rc = (int)floorf(ic/(float)Mc);
            // TODO integer-division

            // column offset in the coarse array
            int cc = ic - rc*Mc;
            // index in the fine array (center element)
            int jf = 2 * rc * Mf + 2*cc;


            // create stencil  // what a waste of precious registers
            // first index: 0=bottom, 1=center, 2=top
            // second index: 0=left, 1=center, 2=right
            float coeffs[3][3] = {0.25f, 0.5f, 0.25f,
                0.5f, 1.0f, 0.5f,
                0.25f, 0.5f, 0.25f};

            int isDirichlet = 0;

            // boundary treatment:
            // simple idea: read in some crap out-of-bounds (rely on zero-padding),
            // and set the corresponding coefficient to zero

            // while we're at it, also figure out if the boundary node
            // has Dirichlet BC, and avoid computations completely
            // (note that the multigrid relies on entries being explicitly
            // set to zero that correspond to Dirichlet boundaries)
            if (rc == 0)
            {
                //bottom
                coeffs[0][0] = coeffs[0][1] = coeffs[0][2] = 0.0f;
                // bottom left node
                if (cc == 0 && n1 == 2) isDirichlet=1;
                // bottom right node
                else if (cc == Mc-1 && n2 == 2) isDirichlet=1;
                // bottom edge
                else if (cc != 0 && cc != Mc-1 && e1 == 2) isDirichlet=1;
            } else if (rc == Mc-1) {

                // top
                coeffs[2][0] = coeffs[2][1] = coeffs[2][2] = 0.0f;
                // top left node
                if (cc == 0 && n4 == 2) isDirichlet=1;
                // top right node
                else if (cc == Mc-1 && n3 == 2) isDirichlet=1;
                // top edge
                else if (cc != 0 && cc != Mc-1 && e3 == 2) isDirichlet=1;
            }

            if (cc == 0)
            {
                //left
                coeffs[0][0] = coeffs[1][0] = coeffs[2][0] = 0.0f;

                // bottom left node
                if (rc == 0 && n1 == 2) isDirichlet=1;
                // top left node
                else if (rc == 0 && n4 == 2) isDirichlet=1;
                // left edge
                else if (rc != 0 && rc != Mc-1 && e4 == 2) isDirichlet=1;
            } else if (cc == Mc-1) {

                // right
                coeffs[0][2] = coeffs[1][2] = coeffs[2][2] = 0.0f;

                // bottom right node
                if (rc == Mc-1 && n2 == 2) isDirichlet=1;
                // top right node
                else if (rc == Mc-1 && n3 == 2) isDirichlet=1;
                // right edge
                else if (rc != 0 && rc != Mc-1 && e2 == 2) isDirichlet=1;
            }


            // all preparations done, so compute if necessary
            if (isDirichlet > 0)
                coarse[ic] = 0.0;
            else
            {
                // compute weighted sum

                // center
                double result = coeffs[1][1] * fine[jf];
                // direct neighbours
                result += coeffs[1][0] * fine[jf-1];
                result += coeffs[1][2] * fine[jf+1];
                result += coeffs[0][1] * fine[jf-Mf];
                result += coeffs[2][1] * fine[jf+Mf];
                // diagonal neighbours
                result += coeffs[0][0] * fine[jf-Mf-1];
                result += coeffs[0][2] * fine[jf-Mf+1];
                result += coeffs[2][0] * fine[jf+Mf-1];
                result += coeffs[2][2] * fine[jf+Mf+1];
                // done
                coarse[ic] = result;
            }

        }
    }
}

extern "C" void cuda_restriction_float(void * coarse, unsigned long size_coarse, void * fine, unsigned long size_fine,
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

    // version directly including Dirichlet treatment
    grid.x = (int)ceil(Nc/(double)block.x);
    honei::cuda::restriction_gpu<<<grid,block>>>(coarse_gpu, fine_gpu, Nc, Nf, Mc, Mf,
            macroBorderMask[0], macroBorderMask[1], macroBorderMask[2], macroBorderMask[3],
            macroBorderMask[4], macroBorderMask[5], macroBorderMask[6], macroBorderMask[7]);
    CUDA_ERROR();

}
