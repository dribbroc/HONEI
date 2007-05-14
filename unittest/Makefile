#Makefile Unittest
.PHONY: all
all: unittest

Main.o	:	Main.cpp BlasTest.hh LiblaTest.hh BasisTest.hh Makefile
		g++ -c Main.cpp -o Main.o

BlasTest.o	:	BlasTest.cpp BlasTest.hh BasisTest.hh cblas.h Makefile
		g++ -c BlasTest.cpp -o BlasTest.o

LiblaTest.o	:	LiblaTest.cpp LiblaTest.hh BasisTest.hh Makefile
		g++ -c LiblaTest.cpp -o LiblaTest.o

unittest	:	Main.o BlasTest.o LiblaTest.o Makefile
		g++ BlasTest.o LiblaTest.o Main.o -o unittest -lg2c -lcblas -lblas -Llib

clean		:	Makefile
		rm *.o
		
