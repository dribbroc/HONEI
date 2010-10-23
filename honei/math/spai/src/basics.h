/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __basics_H
#define __basics_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef MPI
#include "mpi.h"
#define SPAI_Comm MPI_Comm
#else
#define SPAI_Comm void*
#endif

/**********************************************************************/

#ifndef MIN
#  define  MIN(a,b)	((a) < (b) ? (a) : (b))
#endif

void init_SPAI();

void basic_distribution
(SPAI_Comm,
 int,
 int,
 int *,
 int *,
 int *,
 int *,
 int *);

void* mmalloc
(size_t size);

#endif
