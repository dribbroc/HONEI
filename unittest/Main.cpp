#include <iostream>
#include <string>
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
		cout << this->name << endl;
		return true;
	}
	
};

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
BasisTest btest("hallo");
cout << btest.starteTest()<< endl;	
	
	return 0;
	
}





 