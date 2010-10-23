/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __qr_H
#define __qr_H

#include "basics.h"
#include "index_set.h"
#include "spai.h"
#include "spai_error.h"

extern char *Tchar;
extern char *Nchar;
extern char *Uchar;
extern char *Lchar;

#ifdef T3D
extern _fcd Tchar_fcd;
extern _fcd Nchar_fcd;
extern _fcd Uchar_fcd;
extern _fcd Lchar_fcd;
#endif

int qr
(matrix *,
 int,
 int,
 int);

int seek_ptr
(matrix *,
 struct index_set *,
 int,
 int *);

void write_unblocked
(double *v,
 int,
 int,
 int);

void multq
(char *,
 int,
 double *,
 int,
 int);

void fill_zeros
(int,
 int,
 int,
 double *);

#endif
