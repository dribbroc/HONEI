#include "benchmark.hh"
 
using namespace std;
 
// die Prozedur, die in diesem Beispiel gebenchmarkt wird
void zmul(int z, float f1, float f2)
{
    float acc  = 0;
    float temp = 0;
    for (int i = z; i > 0; --i)

    {
        temp = f1*f2;
        acc += temp;
    }

    cout <<  "Ergebnis: " << acc << endl;
}
 
class BenchmarkTest :
    public Benchmark
{
    private:
        int _anzahl;
    public:
        BenchmarkTest(const std::string & type, int anz) :
            Benchmark(type)
        {
            _anzahl = anz;
        }
 
/* es werden 10 mal 2 Zufallszahlen generiert, die _anzahl (50000000) mal in zmul() miteinander multipliziert werden. */
        virtual void run()
        {
            for (int i = 0; i < 10; ++i)  
            {
                float a = rand();
                float b = rand();
                BENCHMARK(zmul(_anzahl, a, b));
            }
            evaluate(2*_anzahl);  //Auswertung generieren.
        }
};
BenchmarkTest TestBench("zmul 50 Mio.", 50000000);
