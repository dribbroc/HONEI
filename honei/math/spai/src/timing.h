
#ifndef __timing_H
#define __timing_H

#include "basics.h"

#define NUMTIMERS 100

extern double timer0[],timer1[];
extern double max_time[],sum_time[];

#define ident_com_server            0

#define ident_handle_get_line       1
#define ident_handle_request_Mline  2
#define ident_handle_put_Mline      3
#define ident_handle_Im_done        4
#define ident_handle_done_signal    5

#define ident_get_line              6
#define ident_request_Mline         7
#define ident_put_Mline             8
#define ident_say_Im_done           9
#define ident_check_done           10

#define ident_block_matrix         11
#define ident_scalar_matrix        12

#define ident_spai                 13
#define ident_bicgstab             14

#define ident_read_mm_matrix       15

void init_timers
();

void start_timer
(int);

void stop_timer
(int);

void report_times
(int,
 char *,
 int,
 SPAI_Comm);

#endif
