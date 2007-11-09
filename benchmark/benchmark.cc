/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.hh>

#include <cstdlib>
#include <utility>
#include <fstream>

// using namespace pg512;
using namespace std;

class Benchmark;

class BenchmarkList
{
    private:
        static std::list<Benchmark *> _tests;

        BenchmarkList()
        {
        }

    public:
        typedef std::list<Benchmark*>::const_iterator Iterator;

        static BenchmarkList * instance()
        {
            static BenchmarkList result;
            return &result;
        }

        void register_test(Benchmark * const test)
        {
            _tests.push_back(test);
        }

        Iterator begin_tests() const
        {
            return _tests.begin();
        }

        Iterator end_tests() const
        {
            return _tests.end();
        }
        
        int bench_count() const
        {
            return _tests.size();
        }
};

std::list<Benchmark *> BenchmarkList::_tests;

Benchmark::Benchmark(const std::string & id) :
    _id(id)
{
    BenchmarkList::instance()->register_test(this);
}

void Benchmark::evaluate()
{
    int x(0), xmin(1), xmax(1);
    double total = 0;
    double min = _benchlist.front();
    double max = _benchlist.front();
    for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
    {
        total += *i;
        ++x;
        if (*i < min)
        {  
            min = *i;
            xmin = x;
        }
        else if (*i > max) 
        {
            max = *i;
            xmax = x;
        } 
    }
    cout << "Function Calls: " << x << endl;
    cout << "Runtime - total: " << total << "sec" << endl;
    cout << "Runtime - lowest: " << min << "sec (" << xmin << ".)" << endl;
    cout << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
    double avg = (total/x);
    cout << "Runtime - average: " << avg << "sec" << endl;
    ofstream ofs("BenchmarkOut.txt", ios_base::out | ios_base::app);
    if (!ofs)
        cout << "Can't write to file!" << endl;
    else
    {
        time_t t;
        time(&t);
        ofs << _id  << " - " << ctime(&t) << endl << endl;
        ofs << "Result:"<< endl;
        ofs << "Function Calls: " << x << endl;
        ofs << "Runtime - total: " << total << "sec" << endl;
        ofs << "Runtime - lowest: " << min << "sec (" << xmin << ".)" << endl;
        ofs << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
        ofs << "Runtime - average: " << avg << "sec" << endl;
        ofs << endl << endl << endl;
    }
}

void Benchmark::evaluate(BenchmarkInfo info)
{
    int x(0), xmin(1), xmax(1);
    double total = 0;
    double min = _benchlist.front();
    double max = _benchlist.front();
    for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
    {
        total += *i;
        ++x;
        if (*i < min)
        {  
            min = *i;
            xmin = x;
        }   
        else if (*i > max)
        { 
            max = *i;
            xmax = x;
        }
    }
    cout << "Function Calls: " << x << endl;
    cout << "Runtime - total: " << total << "sec" << endl;
    cout << "Runtime - lowest: " << min << "sec (" << xmin << ".)" << endl;
    cout << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
    double avg = (total/x);
    cout << "Runtime - average: " << avg << "sec" << endl; 
    string pf = " KMGTPEZY";
    double tp = ((double)(info.load + info.store) / (1024 * 1024)) * x / total;
    int i = 2;
    while (tp > 1024 && i < 8)
    {
        tp /= 1024;
        ++i;
    }
    cout << "Transfer rate: " << tp << pf[i] << "B/s" << endl;
    double f = 0;
    if (total > 0)
        f = ((x/total)*info.flops);
    int j = 0;
    while (f > 1000 && j < 8)
    {
        f /= 1000;
        ++j;
    }
    cout << f << " " << pf[j] << "FLOPS" << endl;
    ofstream ofs("BenchmarkOut.txt", ios_base::out | ios_base::app);
    if (!ofs)
        cout << "Can't write to file!" << endl;
    else
    {
        time_t t;
        time(&t);
        ofs << _id  << " - " << ctime(&t) << endl << endl;
        ofs << "Result:"<< endl;
        ofs << "Function Calls: " << x << endl;
        ofs << "Runtime - total: " << total << "sec" << endl;
        ofs << "Runtime - lowest: " << min << "sec (" << xmin << ".)" << endl;
        ofs << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
        ofs << "Runtime - average: " << avg << "sec" << endl;
        ofs << "Transfer rate: " << tp << pf[i] << "B/s"  << endl;
        ofs << f << " " << pf[j] << "FLOPS" << endl;
        ofs << endl << endl << endl;
    }
}

const std::string Benchmark::id() const
{
    return _id;
}

int main(int argc, char** argv)
{
    int result=EXIT_SUCCESS;
    list<int> runrs;
    cout << "Select Benchmark you'd like to add to runlist:" << endl;
    int count = 1;
    for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_tests()),i_end(BenchmarkList::instance()->end_tests()) ; i != i_end ; ++i, ++count)
    {
        cout << count << ": " << (*i)->id() << endl;
    }
    cout << endl << "Choose: (1), (2), ..., (a)ll, (n)one" << endl;
    string a, tmp;
    int b;
    while (a != "n")
    {
        cout << "Selection: ";
        cin >> a;
        cout << "Added: ";
        if (a == "a")
        {
            cout << "all" << endl;
            a = "n";
            runrs.clear();
            for (int i = 0 ; i < BenchmarkList::instance()->bench_count() ; ++i)
            {
                runrs.push_back(i);
            }
        }
        while ((a.length()>0) && (a != "n"))
        {
            tmp = "";
            b = a[0];
            if ((b < 48) || (b > 57))
                a.erase(0,1);
            else
            {
                while ((b >=  48) && (b <= 65))
                {
                    tmp += a[0];
                    a.erase(0,1);
                    b = a[0];
                } 
                b = atoi(tmp.c_str());
                if ((b > 0) && (b <= BenchmarkList::instance()->bench_count())) 
                {
                    runrs.push_back(b-1);
                    cout << b << " ";
                }
            }
        }
        cout << endl;
        if (a != "n") cout << "Add more? ";
    }
    count = 0;
    if (runrs.size()>0)
    {
        runrs.sort();
        int next = runrs.front();
        runrs.pop_front();
        for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_tests()),i_end(BenchmarkList::instance()->end_tests()) ; i != i_end ; ++i, ++count)
        {
            if (next == count)
            {
                try
                {
                    std::cout << (*i)->id() << ": " << std::endl << std::endl;
                    (*i)->run();
                    std::cout << "'" << (*i)->id() << "' finished successfully!" << endl << endl;
                }
                catch (BenchFailedException & e)
                {
                    std::cout << " FAILED" << std::endl;
                    std::cout << e.what() << std::endl;
                    result = EXIT_FAILURE;
                }
                while ((runrs.size() > 0) && (next == runrs.front()))
                {
                    runrs.pop_front();
                }
                if (runrs.size() > 0)
                {
                    next = runrs.front();
                    runrs.pop_front();
                }
            }
        }
    }
    return result;
}
