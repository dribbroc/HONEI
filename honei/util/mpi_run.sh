#!/bin/bash
# vim: set ft=sh sw=4 sts=4 et :

echo ">>> test ${1} [BACKENDS=${BACKENDS},TYPE=${TYPE}]"
mpirun -np 4 ${@} ${TYPE} ${BACKENDS}
