#include <iostream>
#include <string>
#include "BasisTest.hh"
#include "BlasTest.hh"
using namespace std;


int main(int argc, char** argv) {

/*cout << "Hallo du!\n";
cout << argc << " Argumente. Wert:  " << argv[1];
cout << "\n\n";

cout << "Wieviele Parameter ausgeben?" << endl;

int anzahl;
cin >> anzahl;
cout << anzahl << endl;
	
for (int i = 0; (i < anzahl) && (i < argc); i++) {
	cout << argv[i] << endl; 
} */
BlasTest btest("blub");
//cout << btest.starteTest()<< endl;
btest.starteTest();	
	
	return 0;
	
}





 
