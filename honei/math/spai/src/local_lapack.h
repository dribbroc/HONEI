/* ============================================================================================ *
 *												*
 *  local_lapack.h -- C-Prototypes for the utilised LAPACK subroutines.				*
 *												*
 *  Author: 	Michael Hagemann								*
 *  Created:	2005-10-21 15:02								*
 *												*
 *  Copyright (c) 2005, Michael Hagemann.							*
 *												*
 * -------------------------------------------------------------------------------------------- *
 *												*
 * ============================================================================================ */

/* -------------------------------------------------------------------------------------------- */

#ifndef __LOCAL_LAPACK_H
#define __LOCAL_LAPACK_H


#include "config.h"


#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------------------------- */

	void F77_FUNC(dgetrf, DGETRF) (
		int	*m,
		int	*n,
		double	 a[],
		int	*lda,
		int	 ipiv[],
		int	*info);

	void F77_FUNC(dgetri, DGETRI) (
		int	*n,
		double	 a[],
		int	*lda,
		int	 ipiv[],
		double	 work[],
		int	*lwork,
		int	*info);

	void F77_FUNC(dgeqrf, DGEQRF) (
		int	*m,
		int    	*n,
		double	 a[],
		int    	*lda,
		double	*tau,
		double  *work,
		int    	*lwork,
		int	*info);

	void F77_FUNC(dtrtrs, DTRTRS) (
		char    *uplo,
		char    *trans,
		char    *diag,
		int     *n,
		int     *nrhs,
		double   a[],
		int     *lda,
		double  *b,
		int     *ldb,
		int     *info);

	void F77_FUNC(dormqr, DORMQR) (
		char 	*side,
		char 	*trans,
		int 	*m,
		int 	*n,
		int 	*k,
		double 	 a[],
		int 	*lda,
		double 	*tau,
		double 	 c[],
		int 	*ldc,
		double 	*work,
		int 	*lwork,
		int 	*info);

/* -------------------------------------------------------------------------------------------- */

#ifdef __cplusplus
}
#endif

#endif  /* __LOCAL_LAPACK_H */
