#Makefile Unittest
.PHONY: all
all: unittest

Unittest.o	:	Unittest.cpp BlasTest.hh LiblaTest.hh BasisTest.hh Makefile
		g++ -c Unittest.cpp -o Unittest.o

BlasTest.o	:	BlasTest.cpp BlasTest.hh BasisTest.hh cblas.h Makefile
		g++ -c BlasTest.cpp -o BlasTest.o

LiblaTest.o	:	LiblaTest.cpp LiblaTest.hh BasisTest.hh Makefile
		g++ -c LiblaTest.cpp -o LiblaTest.o

unittest	:	Unittest.o BlasTest.o LiblaTest.o Makefile
		g++ BlasTest.o LiblaTest.o Unittest.o -o unittest -lg2c -lcblas -lblas -Llib

clean		:	Makefile
		rm *.o
		
