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
        static std::list<Benchmark *> _benchs;

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

        void register_bench(Benchmark * const bench)
        {
            _benchs.push_back(bench);
        }

        Iterator begin_benchs() const
        {
            return _benchs.begin();
        }

        Iterator end_benchs() const
        {
            return _benchs.end();
        }
        
        int bench_count() const
        {
            return _benchs.size();
        }

        Iterator erase(int index)
        {
            std::list<Benchmark*>::iterator er(_benchs.begin()), benchs_end(_benchs.end());
            for (int i(1) ; ((i < index) && (er != benchs_end)) ; ++i) 
                ++er;
            if (er != benchs_end)
                return _benchs.erase(er);
            else
                return er;               
        }
};

std::list<Benchmark *> BenchmarkList::_benchs;

Benchmark::Benchmark(const std::string & id) :
    _id(id)
{
    BenchmarkList::instance()->register_bench(this);
}

void Benchmark::calculate()
{
    _x = 0;
    _xmin = 1;
    _xmax = 1;
    _total = 0;
    _min = _benchlist.front();
    _max = _benchlist.front();
    for (list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
    {
        _total += *i;
        ++_x;
        if (*i < _min)
        {  
            _min = *i;
            _xmin = _x;
        }
        else if (*i > _max) 
        {
            _max = *i;
            _xmax = _x;
        } 
    }
    _avg = (_total/_x);
    list<double> temp;
    temp = _benchlist;
    temp.sort();
    list<double>::iterator m = temp.begin();
    for(int i = 0 ; i < ((_x + 1) / 2) ; ++i)
    {
        ++m;
    }
    _median = *m;
}

void Benchmark::calculate(BenchmarkInfo info)
{
    calculate();
    if (_total > 0)
    {
        _tp = ((double)(info.load + info.store) / (1024 * 1024)) * _x / _total;
        _f = ((double)((_x/_total)*info.flops) / 1000000);
    }
    else
    {
        _tp = 0;
        _f = 0;
    }
    if (_median > 0)
    {
        _mediantp = ((double)(info.load + info.store) / (1024 * 1024)) / _median;
        _medianf = ((double)(info.flops) / 1000000) / _median;
    }
    else
    {
        _mediantp = 0;
        _medianf = 0;
    }
}

void Benchmark::evaluate()
{
    calculate();
    cout << "Function Calls: " << _x << endl;
    cout << "Runtime - total:   " << _total << "sec" << endl;
    cout << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << endl;
    cout << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << endl;
    cout << "Runtime - mean:    " << _avg << "sec" << endl;
    cout << "Runtime - median:  " << _median << "sec" << endl;
    ofstream ofs("BenchmarkOut.txt", ios_base::out | ios_base::app);
    if (!ofs)
        cout << "Can't write to file!" << endl;
    else
    {
        time_t t;
        time(&t);
        ofs << _id  << " - " << ctime(&t) << endl << endl;
        ofs << "Result:"<< endl;
        ofs << "Function Calls: " << _x << endl;
        ofs << "Runtime - total: " << _total << "sec" << endl;
        ofs << "Runtime - lowest: " << _min << "sec (" << _xmin << ".)" << endl;
        ofs << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << endl;
        ofs << "Runtime - mean:    " << _avg << "sec" << endl;
        ofs << "Runtime - median:  " << _median << "sec" << endl;
        ofs << endl << endl << endl;
    }
}

void Benchmark::evaluate(BenchmarkInfo info)
{
    calculate(info);
    cout << "Function Calls: " << _x << endl;
    cout << "Runtime - total:   " << _total << "sec" << endl;
    cout << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << endl;
    cout << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << endl;
    cout << "Runtime - mean:    " << _avg << "sec" << endl;
    cout << "Runtime - median:  " << _median << "sec" << endl;
    string pf = " KMGTPEZY";
    int i = 2;
    while (_tp > 1024 && i < 8)
    {
        _tp /= 1024;
        _mediantp /= 1024;
        ++i;
    }
    cout << "Transfer rate (mean):   " << _tp << pf[i] << "B/s" << endl;
    cout << "Transfer rate (median): " << _mediantp << pf[i] << "B/s" << endl;
    int j = 2;
    while (_f > 1000 && j < 8)
    {
        _f /= 1000;
        _medianf /= 1000;
        ++j;
    }
    cout << _f << " " << pf[j] << "FLOPS (mean)" << endl;
    cout << _medianf << " " << pf[j] << "FLOPS (median)" << endl;
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
        ofs << "Function Calls: " << _x << endl;
        ofs << "Runtime - total:   " << _total << "sec" << endl;
        ofs << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << endl;
        ofs << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << endl;
        ofs << "Runtime - mean:    " << _avg << "sec" << endl;
        ofs << "Runtime - median:  " << _median << "sec" << endl;
        ofs << "Transfer rate (mean):   " << _tp << pf[i] << "B/s" << endl;
        ofs << "Transfer rate (median): " << _mediantp << pf[i] << "B/s" << endl;
        ofs << _f << " " << pf[j] << "FLOPS (mean)" << endl;
        ofs << _medianf << " " << pf[j] << "FLOPS (median)" << endl;
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
            ofs << std::setw(13) << _medianf << "\t" << std::setw(11) << _mediantp << "\t" << std::setw(11) << _min << "\t" << std::setw(11) << _max << "\t" << std::setw(12) << _avg << "\t" << std::setw(14) << _median << "\t" << std::setw(11) << _f << "\t" << std::setw(9) << _tp;
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

void Benchmark::register_tag(std::string tag_name)
{
    _tag_name = tag_name;
}

std::string Benchmark::get_tag_name()
{
    return _tag_name;
}

int main(int argc, char** argv)
{
    int result=EXIT_SUCCESS;
    list<int> runrs;
    bool sse(true);
    bool cell(true);
    bool mc(true);
    bool sc(true);
    if (argc > 1)
    {
        sse = false;
        cell = false;
        mc = false;
        sc = false;
        for(int i(1) ; i < argc ; ++i)
        {
            if (honei::stringify(argv[i]) == "sse")
            {
                sse = true;
            }
            if (honei::stringify(argv[i]) == "cell")
            {
                cell = true;
            }
            if (honei::stringify(argv[i]) == "mc")
            {
                mc = true;
            }
            if (honei::stringify(argv[i]) == "cpu")
            {
                sse = true;
                mc = true;
                sc = true;
            }
            if (honei::stringify(argv[i]) == "sc")
            {
                sc = true;
            }
        }
    }
    cout << "Select Benchmark you'd like to add to runlist:" << endl;
    int count = 1;
    for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_benchs()),i_end(BenchmarkList::instance()->end_benchs()) ; i != i_end ; )
    {
        if (sse && (((*i)->get_tag_name() == "sse") || ((*i)->get_tag_name() == "mc-sse")))
        {
            cout << count << ": " << (*i)->id() << endl;
            ++count;
            ++i;
            continue;
        }
        if (cell && ((*i)->get_tag_name() == "cell"))
        {
            cout << count << ": " << (*i)->id() << endl;
            ++count;
            ++i;
            continue;
        }
        if (mc && ((*i)->get_tag_name() == "mc"))
        {
            cout << count << ": " << (*i)->id() << endl;
            ++count;
            ++i;
            continue;
        }
        if (sc && ((*i)->get_tag_name() == "cpu"))
        {
            cout << count << ": " << (*i)->id() << endl;
            ++count;
            ++i;
            continue;
        }
        i = BenchmarkList::instance()->erase(count);

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
        for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_benchs()),i_end(BenchmarkList::instance()->end_benchs()) ; i != i_end ; ++i, ++count)
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
