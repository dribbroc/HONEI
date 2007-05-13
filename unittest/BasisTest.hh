#ifndef BASISTEST_HH
#define BASISTEST_HH 1
#include <string>
#include <iostream>
using namespace std;	
class BasisTest  {
	string name;
	string beschreibung;

	public:
	
	BasisTest(string n) {
		this->name         = n;
		this->beschreibung = "Ein Basistest, gibt immer true zurück!";	
	}


	bool starteTest() {
		std::cout << this->name << endl;
		return true;
	}
	
	
	
};
#endif
