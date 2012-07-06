#!/bin/bash
# @ job_name = arm
# @ initialdir = .
# @ output = arm_%j.out
# @ error = arm_%j.err
# @ total_tasks = 192
# @ cpus_per_task = 1
# @ wall_clock_limit = 00:30:00
#mpiexec -n 64 sh -c 'echo 1000 >> /proc/self/oom_score_adj; ./honei-mpi-ring 2500 2500 25 config'

srun honei-mpi-ring 2000 2000 25 config
