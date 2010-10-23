
#include <stdio.h>

#include "mex.h"
#include "spai.h"
#include "read_mm_matrix.h"


/* --------------------------------------------------------------------- */

matrix	*matlab_to_matrix   (mxArray	 *Aml);
mxArray	*matrix_to_matlab   (matrix	 *M);
void	 matlab_to_mm_data  (mxArray  	 *Aml,
			     mm_data 	**rows_ptr,
			     int 	 *nnz_rows_ptr,
			     mm_data 	**cols_ptr,
			     int 	 *nnz_cols_ptr,
			     int 	 *N_ptr);


/* --------------------------------------------------------------------- */

mxArray *spai_full (mxArray *Aml,
		    double   ep,
		    int      ns,
		    int      mb,
		    int      mn,
		    int      bs,
		    int      vb,
		    int      sc,
		    int      ld,
		    int      ud,
		    double   tau)
{
    matrix  *M, *A;
    mxArray *Mml;
    int      cs = 0;
    int      status;

    if (! mxIsSparse(Aml))
	mexErrMsgTxt ("The array is not sparse.");

    if (! mxIsNumeric(Aml))
	mexErrMsgTxt ("The array is not numeric.");

    if (mxIsComplex(Aml))
	mexErrMsgTxt ("The array is not real.");

    if (mxGetN(Aml) != mxGetM(Aml))
    {
	mexErrMsgTxt ("The array is not square.");
    }

    A = matlab_to_matrix (Aml);

    status = bspai (A, &M, stderr,
		    ep, ns, mb, mn, bs, cs, vb,
		    sc, ld, ud, tau);

    Mml = matrix_to_matlab (M);

    sp_free_matrix (A);
    sp_free_matrix (M);

    return (Mml);
}


/* --------------------------------------------------------------------- */

matrix *matlab_to_matrix (mxArray *Aml)
{
    matrix  *A;
    mm_data *rows, *cols;
    int      N, nnz_rows, nnz_cols, mnl;

    matlab_to_mm_data (Aml,
		       &rows,
		       &nnz_rows,
		       &cols,
		       &nnz_cols,
		       &N);

    mnl = N;

    A = mm_to_matrix (0,
		      0,
		      N,
		      1,
		      mnl,
		      rows, nnz_rows,
		      cols, nnz_cols,
		      NULL);

    free (rows);
    free (cols);

    return (A);
}


/* --------------------------------------------------------------------- */

void matlab_to_mm_data (mxArray  *Aml,
			mm_data **rows_ptr,
			int 	 *nnz_rows_ptr,
			mm_data **cols_ptr,
			int 	 *nnz_cols_ptr,
			int 	 *N_ptr)
{
    mm_data *rows, *cols;
    int      N, M, nnz_rows, nnz_cols;
    int     *ir, *jc;
    double  *pr;
    int      len, j, k;

    M  = N = mxGetN (Aml);
    pr = mxGetPr (Aml);
    ir = mxGetIr (Aml);
    jc = mxGetJc (Aml);
    nnz_rows = nnz_cols = jc[N];

    *N_ptr 	  = N;
    *nnz_rows_ptr = nnz_rows;
    *nnz_cols_ptr = nnz_cols;

    rows = (mm_data *) malloc (nnz_rows * sizeof(mm_data));
    cols = (mm_data *) malloc (nnz_cols * sizeof(mm_data));

    *rows_ptr = rows;
    *cols_ptr = cols;

    len = 0;

    for (j = 0; j  <  N; j++)
    {
	for (k = jc[j]; k  <  jc[j+1]; k++)
	{
	    rows[len].i   = ir[k];
	    rows[len].ib  = 0;
	    rows[len].j   = j;
	    rows[len].jb  = 0;
	    rows[len].val = pr[k];

	    cols[len].i   = ir[k];
	    cols[len].ib  = 0;
	    cols[len].j   = j;
	    cols[len].jb  = 0;
	    cols[len].val = pr[k];

	    len++;
	}
    }
}


/* --------------------------------------------------------------------- */

mxArray *matrix_to_matlab (matrix *A)
{
    mxArray  *Mml_full;
    double   *pr;
    int       N, nnz, j, i, k, col, row, bs, bs2, ib, jb, index;

    N 	     = A->n;
    bs 	     = 1;
    bs2      = bs * bs;
    nnz      = count_nonzeros(A) * bs2;
    Mml_full = mxCreateDoubleMatrix (nnz, 3, 0);
    pr 	     = mxGetPr (Mml_full);

    for (j = 0, k = 0;  j < A->mnls[0];  j++)
    {
	for (i = 0;  i < A->lines->len[j];  i++)
	{
	    index = A->lines->ptrs[j][i];

	    if (A->transposed == 0)
	    {
		row = index;
		col = j + A->start_indices[0];
	    }
	    else
	    {
		col = index;
		row = j + A->start_indices[0];
	    }

	    row = 1 + bs * row;
	    col = 1 + bs * col;

	    for (ib = 0;  ib < bs;  ib++)
	    {
		for (jb = 0;  jb < bs;  jb++)
		{
		    pr[k] 	  = row+ib;
		    pr[k+nnz]     = col+jb;
		    pr[k+nnz+nnz] = A->lines->A[j][bs2*i + jb*bs + ib];
		    k++;
		}
	    }
	}
    }

    return (Mml_full);
}


/* --------------------------------------------------------------------- */

/* Gateway function. */

void mexFunction (int		 nlhs,
		  mxArray	*plhs[],
                  int		 nrhs,
		  const mxArray	*prhs[])
{
    double 	 ep;
    int 	 ns;
    int 	 mb;
    int 	 mn;
    int 	 bs;
    int 	 vb;
    int 	 sc;
    int 	 ld;
    int 	 ud;
    double 	 ta;

    mxArray *Aml;
    mxArray *Mml;

    if (nrhs == 0)
    {
	mexErrMsgTxt ("spai_full requires at least one input arg.");
    }

    Aml = (mxArray *) prhs[0];
    ep  = *mxGetPr (prhs[1]);
    ns  = *mxGetPr (prhs[2]);
    mb  = *mxGetPr (prhs[3]);
    mn  = *mxGetPr (prhs[4]);
    bs  = *mxGetPr (prhs[5]);
    vb  = *mxGetPr (prhs[6]);
    sc  = *mxGetPr (prhs[7]);
    ld  = *mxGetPr (prhs[8]);
    ud  = *mxGetPr (prhs[9]);
    ta  = *mxGetPr (prhs[10]);

    if (vb)
    {
	printf ("SPAI parameters:\n\n");
	printf ("   ep = %.2le\n", ep);
	printf ("   ns = %d\n", ns);
	printf ("   mn = %d\n", mn);
	printf ("   bs = %d\n", bs);
	printf ("   mb = %d\n", mb);
	printf ("   vb = %d\n", vb);
	printf ("   sc = %d\n", sc);
	printf ("   ld = %d\n", ld);
	printf ("   ud = %d\n", ud);
	printf ("   ta = %.2le\n", ta);
    }

    plhs[0] = spai_full(Aml, ep, ns, mb, mn, bs, vb, sc, ld, ud, ta);

    return;
}
