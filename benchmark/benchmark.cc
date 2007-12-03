/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.hh>

#include <cstdlib>
#include <utility>
#include <iomanip>
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

void Benchmark::calculate()
{
    x = 0;
    xmin = 1;
    xmax = 1;
    total = 0;
    min = _benchlist.front();
    max = _benchlist.front();
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
    avg = (total/x);
    list<double> temp;
    temp = _benchlist;
    temp.sort();
    list<double>::iterator m = temp.begin();
    for(int i = 0 ; i < ((x + 1) / 2) ; ++i)
    {
        ++m;
    }
    median = *m;
}

void Benchmark::calculate(BenchmarkInfo info)
{
    calculate();
    if (total > 0)
    {
        tp = ((double)(info.load + info.store) / (1024 * 1024)) * x / total;
        f = ((double)((x/total)*info.flops) / 1000000);
    }
    else
    {
        tp = 0;
        f = 0;
    }
    if (median > 0)
    {
        mediantp = ((double)(info.load + info.store) / (1024 * 1024)) / median;
        medianf = ((double)(info.flops) / 1000000) / median;
    }
    else
    {
        mediantp = 0;
        medianf = 0;
    }
}

void Benchmark::evaluate()
{
    calculate();
    cout << "Function Calls: " << x << endl;
    cout << "Runtime - total:   " << total << "sec" << endl;
    cout << "Runtime - lowest:  " << min << "sec (" << xmin << ".)" << endl;
    cout << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
    cout << "Runtime - mean:    " << avg << "sec" << endl;
    cout << "Runtime - median:  " << median << "sec" << endl;
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
        ofs << "Runtime - mean:    " << avg << "sec" << endl;
        ofs << "Runtime - median:  " << median << "sec" << endl;
        ofs << endl << endl << endl;
    }
}

void Benchmark::evaluate(BenchmarkInfo info)
{
    calculate(info);
    cout << "Function Calls: " << x << endl;
    cout << "Runtime - total:   " << total << "sec" << endl;
    cout << "Runtime - lowest:  " << min << "sec (" << xmin << ".)" << endl;
    cout << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
    cout << "Runtime - mean:    " << avg << "sec" << endl;
    cout << "Runtime - median:  " << median << "sec" << endl;
    string pf = " KMGTPEZY";
    int i = 2;
    while (tp > 1024 && i < 8)
    {
        tp /= 1024;
        mediantp /= 1024;
        ++i;
    }
    cout << "Transfer rate (mean):   " << tp << pf[i] << "B/s" << endl;
    cout << "Transfer rate (median): " << mediantp << pf[i] << "B/s" << endl;
    int j = 2;
    while (f > 1000 && j < 8)
    {
        f /= 1000;
        medianf /= 1000;
        ++j;
    }
    cout << f << " " << pf[j] << "FLOPS (mean)" << endl;
    cout << medianf << " " << pf[j] << "FLOPS (median)" << endl;
    if (info.scale != 1)
    {
        cout << "Dense version calculates " << info.scale << " times more flops.\n(Depends on used elements of Sparse/Banded Operands.)" << endl;
        if (info.scaleinfo != " ")
            cout << info.scaleinfo << endl;
    }
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
        ofs << "Runtime - total:   " << total << "sec" << endl;
        ofs << "Runtime - lowest:  " << min << "sec (" << xmin << ".)" << endl;
        ofs << "Runtime - highest: " << max << "sec (" << xmax << ".)" << endl;
        ofs << "Runtime - mean:    " << avg << "sec" << endl;
        ofs << "Runtime - median:  " << median << "sec" << endl;
        ofs << "Transfer rate (mean):   " << tp << pf[i] << "B/s" << endl;
        ofs << "Transfer rate (median): " << mediantp << pf[i] << "B/s" << endl;
        ofs << f << " " << pf[j] << "FLOPS (mean)" << endl;
        ofs << medianf << " " << pf[j] << "FLOPS (median)" << endl;
        ofs << endl << endl << endl;
    }
}

void Benchmark::evaluate_to_plotfile(std::list<BenchmarkInfo> info, std::list<int> cores, int count)
{
    std::list<double> BenchlistCopy;
    BenchlistCopy = _benchlist;
    list<double>::iterator blc = BenchlistCopy.begin();
    list<BenchmarkInfo>::iterator j = info.begin();
    size_t pos;
    time_t t;
    time(&t);
    std::string filename(std::string("PlotOut ") + ctime(&t));
    while ((pos=filename.find(" "))!=-1) filename.replace(pos, 1, "_");
    while ((pos=filename.find(":"))!=-1) filename.replace(pos, 1, "_");
    filename.replace(32, 1, ".txt");
    ofstream ofs(filename.c_str(), ios_base::out);
    if (!ofs)
        cout << "Can't write to file!" << endl;
    else
    {
        ofs.setf(ios::left, ios::adjustfield);
        ofs << "#" << _id << "\n";
        ofs << "#Cores\t";
        if (!((j->size).empty()))
        {
            int counter(1);
            for(list<unsigned long>::iterator si((j->size).begin()); si != (j->size).end() ; ++si, ++counter)
            {
                ofs << counter <<".operand size\t";
            }
        }
        ofs << "median MFLOPS\tmedian MB/s\tmin runtime\tmax runtime\tmean runtime\tmedian runtime\tmean MFLOPS\tmean MB/s\tlist of all runtimes";
        for(list<int>::iterator i = cores.begin() ; i != cores.end() ; ++i, ++j)
        {
            ofs << "\n" << std::setw(6) << *i << "\t";
            if (!((j->size).empty()))
                for(list<unsigned long>::iterator si((j->size).begin()) ; si != (j->size).end() ; ++si)
                {
                    ofs << std::setw(14) << *si << "\t";
                } 
            _benchlist.clear();
            for(int k(0) ; k < count ;  ++k)
            {
                _benchlist.push_back(*blc);
                ++blc;
            }
            calculate(*j);
            ofs << std::setw(13) << medianf << "\t" << std::setw(11) << mediantp << "\t" << std::setw(11) << min << "\t" << std::setw(11) << max << "\t" << std::setw(12) << avg << "\t" << std::setw(14) << median << "\t" << std::setw(11) << f << "\t" << std::setw(9) << tp;
            for (list<double>::iterator bl = _benchlist.begin() ; bl != _benchlist.end() ; ++bl)
            {
               ofs << "\t" << std::setw(8) << *bl;
            }
        }
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
                    std::cout << (*i)->id() << ": " << std::endl;
                    (*i)->run();
                    std::cout << "'" << (*i)->id() << "' finished successfully!" << std::endl << std::endl << std::endl;
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
