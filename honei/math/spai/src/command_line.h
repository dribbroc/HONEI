/*
   SPAI MPI/C version @Copyright 1996,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __command_line_H
#define __command_line_H

#include <string.h>
#include "basics.h"

#ifndef MPI
#define MPI_Comm void*
#endif

void parameters
(int,
 char **,
 char **,
 char **,
 double *,
 int *,
 int *,
 int *,
 int *,
 char **,
 int *,
 double *,
 int *,
 int *,
 int *,
 int *,
 int *,
 int *,
 int *,
 int *,
 int *,
 double *,
 MPI_Comm);

int match_arg(char *string);

#endif

