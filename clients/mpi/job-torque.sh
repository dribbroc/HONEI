#!/bin/bash

# Output files will be named "job-name.{o,e}[0-9]*"
#PBS -N arm
# Get two full nodes, one process per node
# Please use "procs=N" if you do not care about exclusive node acces
# #PBS -l procs=8
#PBS -l nodes=64:ppn=2
# #PBS -l nodes=8
# Queue name (see "qstat -q" for a list)
#PBS -q batch
# Expected execution time (optional; helps the scheduler improve cluster utilization)
#PBS -l walltime=00:30:00
# Do not re-run the job if it fails
#PBS -r n

# Be nice to the system if we try to allocate more physical memory than it is
# available (otherwise the system could partially crash).
echo 1000 >> /proc/self/oom_score_adj

# Make the given program span all the allocated nodes for the job.
#
# We also need to make mpiexec children behave nicely (w.r.t. OOM).
cd /home/dribbrock/honei/trunk/clients/mpi/
#cat $PBS_NODEFILE
#cat $PBS_NODEFILE > nodes.txt
#cat nodes.txt | awk '{print $0 ":2"}' > hostfile
#echo "================="

#mpiexec -np 2 sh -c 'echo 1000 >> /proc/self/oom_score_adj; echo "I am process $$ (son of $PPID) at `hostname` oom=`cat /proc/self/oom_score`"'
#mpiexec -n 8 -f hostfile sh -c 'echo "hostname is `hostname`"'
#mpiexec -n 8 sh -c 'echo "hostname is `hostname`"'
mpiexec -n 128 sh -c 'echo 1000 >> /proc/self/oom_score_adj; ./honei-mpi-ring 2000 2000 25 config'


# 96 nodes
# dualcore
# wie exklusiv usage? -> tar ball von dom .sh script
# 800/850 MB pro knoten -> 400 pro prozess
# fortran flags: specfem->flags.guess
