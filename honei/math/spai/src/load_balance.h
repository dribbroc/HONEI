/*
   SPAI Version 3.0 @Copyright 1999,  All Rights Reserved
   Stephen Barnard
*/

#ifndef __load_balance_H
#define __load_balance_H

#include "basics.h"
#include "vector.h"
#include "index_set.h"
#include "matrix.h"
#include "spai.h"
#include "com_server.h"

extern int next_line;

int grab_Mline
(matrix *,
 matrix *,
 SPAI_Comm);

#endif
