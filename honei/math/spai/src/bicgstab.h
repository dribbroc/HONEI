/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __bicgstab_H
#define __bicgstab_H

#include "basics.h"
#include "vector.h"
#include "index_set.h"
#include "matrix.h"
#include "timing.h"

#include "debug.h"

int bicgstab_R
(void Av_proc(),
 matrix *,
 matrix *,
 vector *,
 vector *,
 int,
 double,
 int);

int bicgstab_L
(void Av_proc(),
 matrix *,
 matrix *,
 vector *,
 vector *,
 int,
 double,
 int);

#endif

