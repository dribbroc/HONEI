/* vim: set sw=4 sts=4 et foldmethod=syntax : */
#include <iostream>
#include <string>
#include "basetest.hh"
#include "blastest.hh"
#include "liblatest.hh"
using namespace std;


int main(int argc, char** argv) 
{

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
//cout << btest.start_test()<< endl;
btest.start_test();
BlasTest btest2();	
LiblaTest latest();
	
	return 0;
	
}





 
