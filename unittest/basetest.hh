/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#ifndef BASETEST_HH
#define BASETEST_HH 1
#include <string>
#include <iostream>
using namespace std;	
class BaseTest  
{
	string name;
	string beschreibung;

	public:
	
	BaseTest(string n) 
    {
		this->name         = n;
		this->beschreibung = "Ein Basistest, gibt immer true zurück!";	
	}
    /// overwrite value at position, returns the changed vector
    //virtual Vector * vector_rw(Vector vector,int position, double value);
     
	 virtual bool start_test()=0;  // Rein virtuell.

	
	
	
};
#endif
