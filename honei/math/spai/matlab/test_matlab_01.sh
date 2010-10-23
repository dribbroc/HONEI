#!/bin/sh

rm -f test_matlab_01.check
rm -f test_matlab_01.checkraw

matlab -nojvm -nosplash > test_matlab_01.checkraw 2>&1  << END
addpath('./')
load m1.dat
A   = spconvert(m1);
M   = spai(A,0.4,10);
AM  = A*M;
s   = size(A);
rhs = A*ones(s(1),1);
[x,flag,relres,iter] = bicgstab(AM,rhs,1e-9,100)
END

# Remove first 9 lines with Matlab version
cat test_matlab_01.checkraw | sed "1,9cmatlab" > test_matlab_01.check

diff  test_matlab_01.check  test_matlab_01.out

if test $? -ne 0
then
  echo check difference in output!!!
fi
