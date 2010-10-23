#!/bin/sh

rm -f test_spai_01.check

# Run tests, but exclude lines with time and version
./spai ../data/m1.mm ../data/m1_rhs.mm -sc 0 | grep -v "\+ SPAI" >   test_spai_01.check
./spai ../data/m1.mm ../data/m1_rhs.mm -sc 1 | grep -v "\+ SPAI" >>  test_spai_01.check
./spai ../data/m1.mm ../data/m1_rhs.mm -sc 2 | grep -v "\+ SPAI" >>  test_spai_01.check

diff  test_spai_01.check  test_spai_01.out

if test $? -ne 0
then
  echo check difference in output!!!
fi
