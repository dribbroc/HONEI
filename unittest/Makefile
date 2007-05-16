#Makefile Unittest
.PHONY: all
all: unittest

unittest.o	:	unittest.cc blastest.hh liblatest.hh basetest.hh Makefile
		g++ -c unittest.cc -o unittest.o

blastest.o	:	blastest.cc blastest.hh basetest.hh Makefile
		g++ -c blastest.cc -o blastest.o

liblatest.o	:	liblatest.cc liblatest.hh basetest.hh Makefile
		g++ -c liblatest.cc -o liblatest.o

testfactory.o	:	testfactory.cc testfactory.hh Makefile
		g++ -c testfactory.cc -o testfactory.o

testdispatcher.o	:	testdispatcher.cc testdispatcher.hh tests.hh Makefile
		g++ -c testdispatcher.cc -o testdispatcher.o

unittest	:	unittest.o blastest.o liblatest.o testfactory.o testdispatcher.o Makefile
		g++ blastest.o liblatest.o unittest.o testfactory.o testdispatcher.o -o unittest -lg2c -lcblas

clean		:	Makefile
		rm -f *.o
		
