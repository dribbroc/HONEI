/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <benchmark/benchmark.hh>

#include <cstdlib>
#include <utility>
#include <iomanip>
#include <fstream>

// using namespace pg512;

class Benchmark;

class BenchmarkList
{
    private:
        static std::list<Benchmark *> _benchs;

        BenchmarkList()
        {
        }

    public:
        typedef std::list<Benchmark*>::iterator Iterator;

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

        Iterator erase(Iterator index)
        {
            return _benchs.erase(index);             
        }
};

std::list<Benchmark *> BenchmarkList::_benchs;

Benchmark::Benchmark(const std::string & id) :
    _id(id),
    _plots(false)
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
    for (std::list<double>::iterator i = _benchlist.begin() ; i != _benchlist.end() ; ++i)
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
    std::list<double> temp;
    temp = _benchlist;
    temp.sort();
    std::list<double>::iterator m = temp.begin();
    if (temp.size() > 1)
    {
        for(int i = 0 ; i < ((_x + 1) / 2) ; ++i)
        {
            ++m;
        }
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
        _mintp = ((double)(info.load + info.store) / (1024 * 1024)) / _min;
        _medianf = ((double)(info.flops) / 1000000) / _median;
        _minf = ((double)(info.flops) / 1000000) / _min;
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
    std::cout << "Function Calls: " << _x << std::endl;
    std::cout << "Runtime - total:   " << _total << "sec" << std::endl;
    std::cout << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << std::endl;
    std::cout << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << std::endl;
    std::cout << "Runtime - mean:    " << _avg << "sec" << std::endl;
    std::cout << "Runtime - median:  " << _median << "sec" << std::endl;
    /// \todo Reenable file output if needed
    /*std::ofstream ofs("BenchmarkOut.txt", std::ios_base::out | std::ios_base::app);
    if (!ofs)
        std::cout << "Can't write to file!" << std::endl;
    else
    {
        time_t t;
        time(&t);
        ofs << _id  << " - " << ctime(&t) << std::endl << std::endl;
        ofs << "Result:"<< std::endl;
        ofs << "Function Calls: " << _x << std::endl;
        ofs << "Runtime - total: " << _total << "sec" << std::endl;
        ofs << "Runtime - lowest: " << _min << "sec (" << _xmin << ".)" << std::endl;
        ofs << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << std::endl;
        ofs << "Runtime - mean:    " << _avg << "sec" << std::endl;
        ofs << "Runtime - median:  " << _median << "sec" << std::endl;
        ofs << std::endl << std::endl << std::endl;
    }*/
}

void Benchmark::evaluate(BenchmarkInfo info)
{
    calculate(info);
    std::cout << "Function Calls: " << _x << std::endl;
    std::cout << "Runtime - total:   " << _total << "sec" << std::endl;
    std::cout << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << std::endl;
    std::cout << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << std::endl;
    std::cout << "Runtime - mean:    " << _avg << "sec" << std::endl;
    std::cout << "Runtime - median:  " << _median << "sec" << std::endl;
    std::string pf = " KMGTPEZY";
    int i = 2;
    while (_tp > 1024 && i < 8)
    {
        _tp /= 1024;
        _mediantp /= 1024;
        _mintp /= 1024;
        ++i;
    }
    std::cout << "Transfer rate (mean):   " << _tp << " " << pf[i] << "B/s" << std::endl;
    std::cout << "Transfer rate (median): " << _mediantp << " " << pf[i] << "B/s" << std::endl;
    std::cout << "Transfer rate (peak):   " << _mintp << " " << pf[i] << "B/s" << std::endl;
    int j = 2;
    while (_f > 1000 && j < 8)
    {
        _f /= 1000;
        _medianf /= 1000;
        _minf /= 1000;
        ++j;
    }
    std::cout << _f << " " << pf[j] << "FLOPS (mean)" << std::endl;
    std::cout << _medianf << " " << pf[j] << "FLOPS (median)" << std::endl;
    std::cout << _minf << " " << pf[j] << "FLOPS (peak)" << std::endl;
    if (info.scale != 1)
    {
        std::cout << "Dense version calculates " << info.scale << " times more flops.\n(Depends on used elements of Sparse/Banded Operands.)" << std::endl;
        if (info.scaleinfo != " ")
            std::cout << info.scaleinfo << std::endl;
    }
    /// \todo Reenable file output if needed
    /*std::ofstream ofs("BenchmarkOut.txt", std::ios_base::out | std::ios_base::app);
    if (!ofs)
        std::cout << "Can't write to file!" << std::endl;
    else
    {
        time_t t;
        time(&t);
        ofs << _id  << " - " << ctime(&t) << std::endl << std::endl;
        ofs << "Result:"<< std::endl;
        ofs << "Function Calls: " << _x << std::endl;
        ofs << "Runtime - total:   " << _total << "sec" << std::endl;
        ofs << "Runtime - lowest:  " << _min << "sec (" << _xmin << ".)" << std::endl;
        ofs << "Runtime - highest: " << _max << "sec (" << _xmax << ".)" << std::endl;
        ofs << "Runtime - mean:    " << _avg << "sec" << std::endl;
        ofs << "Runtime - median:  " << _median << "sec" << std::endl;
        ofs << "Transfer rate (mean):   " << _tp << " " << pf[i] << "B/s" << std::endl;
        ofs << "Transfer rate (median): " << _mediantp << " " << pf[i] << "B/s" << std::endl;
        ofs << "Transfer rate (peak):   " << _mintp << " " << pf[i] << "B/s" << std::endl;
        ofs << _f << " " << pf[j] << "FLOPS (mean)" << std::endl;
        ofs << _medianf << " " << pf[j] << "FLOPS (median)" << std::endl;
        ofs << _minf << " " << pf[j] << "FLOPS (peak)" << std::endl;
        ofs << std::endl << std::endl << std::endl;
    }*/
}

void Benchmark::evaluate_to_plotfile(std::list<BenchmarkInfo> info, std::list<std::string> cores, int count)
{
    std::list<double> BenchlistCopy;
    BenchlistCopy = _benchlist;
    std::list<double>::iterator blc = BenchlistCopy.begin();
    std::list<BenchmarkInfo>::iterator j = info.begin();
    int sizes(0);
    int nr, cnr(1);
    std::string temp;
    time_t t;
    time(&t);
    if ((cores.size() > 1) && (*(cores.begin()) == *(++cores.begin())) && (*(cores.begin()) != *(--cores.end())))
    {
        temp = *cores.begin();
        for (std::list<std::string>::iterator i = cores.begin() ; temp != *(--cores.end()) ; ++i)
        {
            if (temp != *i)
            {
                ++cnr;
                temp = *i;
            }
        }
    }
    nr = cnr;
    std::ifstream ifs("BenchmarkPlotData", std::ifstream::in);
    if (ifs.is_open())
    {
        ifs.ignore(1, '#');
        ifs >> nr;
        nr += cnr;
    }
    ifs.close();
    if (nr == cnr)
    {
        std::ofstream ofs1("BenchmarkPlotData", std::ios_base::out);
            ofs1 << "#" << nr << "\n";
        ofs1.close();
    }
    else
    {
        std::ofstream ofs1("BenchmarkPlotData", std::ios_base::in);
            ofs1 << "#" << nr;
        ofs1.close();
    }
    std::ofstream ofs3("RecentPlots.tex", std::ios_base::app);
    ofs3 << "\\section{" << _id << "}\n";
    ofs3 << "\t\\begin{center}\n";
    ofs3 << "\t\\begin{longtable}{c c c c c}\n";
    ofs3 << "\t\tX & operand size & median runtime & median MFLOPS & median transferrate\\\\ [0.5ex]\n";
    ofs3 << "\t\t\\hline\n";
    ofs3 << "\t\t\\endfirsthead\n";
    ofs3 << "\t\tX & operand size & median runtime & median MFLOPS & median transferrate\\\\ [0.5ex]\n";
    ofs3 << "\t\t\\hline\n";
    ofs3 << "\t\t\\endhead\n";
    std::ofstream ofs("BenchmarkPlotData", std::ios_base::out | std::ios_base::app);
    if (!ofs)
        std::cout << "Can't write to file!" << std::endl;
    else
    {
        ofs.setf(std::ios::left, std::ios::adjustfield);
        ofs << "#" << _id << " --- " << ctime(&t);
        ofs << "#Cores\t";
        if (!((j->size).empty()))
        {
            sizes = 1;
            for(std::list<unsigned long>::iterator si((j->size).begin()); si != (j->size).end() ; ++si, ++sizes)
            {
                ofs << sizes <<".operand size\t";
            }
        }
        ofs << "median MFLOPS\tmedian MB/s\tmin runtime\tmax runtime\tmean runtime\tmedian runtime\tmean MFLOPS\tmean MB/s\tlist of all runtimes";
        if (cnr > 1)
            temp = *cores.begin();
        for(std::list<std::string>::iterator i = cores.begin() ; i != cores.end() ; ++i)
        {
            if ((cnr > 1) && (temp != *i))
            {
                ofs << "\n\n";
                temp = *i;
            }
            if (j != info.end())
            {
                ofs << "\n" << std::setw(6) << *i << "\t";
                if (!((j->size).empty()))
                    for(std::list<unsigned long>::iterator si((j->size).begin()) ; si != (j->size).end() ; ++si)
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
                ofs3 << "\t\t" << *i << " & " << (j->size).front() << " & " << _median << " & " << _medianf << " & " << _mediantp << " \\\\\n";
                ofs << std::setw(13) << _medianf << "\t" << std::setw(11) << _mediantp << "\t" << std::setw(11) << _min << "\t" << std::setw(11) << _max << "\t" << std::setw(12) << _avg << "\t" << std::setw(14) << _median << "\t" << std::setw(11) << _f << "\t" << std::setw(9) << _tp;
                for (std::list<double>::iterator bl = _benchlist.begin() ; bl != _benchlist.end() ; ++bl)
                {
                    ofs << "\t" << std::setw(8) << *bl;
                }
                ++j;
            }
            else
            {
                ofs3 << "% [PDH] paste " << *i << " data here.\n";
                ofs3 << "\t\t" << *i << " & " << 0 << " & " << 0 << " & " << 0 << " & " << 0 << " \\\\\n";
                ofs << "\n# [PDH] paste " << *i << " data here.\n";
                ofs << std::setw(13) << 0 << "\t" << std::setw(11) << 0 << "\t" << std::setw(11) << 0 << "\t" << std::setw(11) << 0 << "\t" << std::setw(12) << 0 << "\t" << std::setw(14) << 0 << "\t" << std::setw(11) << 0 << "\t" << std::setw(9) << 0;
            }
        }
        ofs << "\n\n\n";
    }
    ofs.close();
    --sizes;
    size_t pos;
    bool plotsvx(not(*(info.begin()->size.begin()) == *((--info.end())->size.begin()))), plotcvx(not(*cores.begin() == *(--cores.end())));
    std::string filename(std::string("PlotOut_") + ctime(&t));
    while ((pos=filename.find(" "))!=-1) filename.replace(pos, 1, "_");
    while ((pos=filename.find(":"))!=-1) filename.replace(pos, 1, "_");
    std::string eps1name(filename), eps2name(filename), eps3name(filename);
    filename.replace(32, 1, ".plt");
    eps1name.replace(32, 1, ".eps");
    eps2name.replace(32, 1, "_2.eps");
    eps3name.replace(32, 1, "_3.eps");
    std::ofstream ofs2(filename.c_str(), std::ios_base::out);
    if (plotsvx && not plotcvx)
    {
        ofs2 << "set terminal postscript eps color\nset key below\nset yrange [0:]\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"time in sec.\"\nset output \"" << eps1name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 7+sizes << " t \"median runtime\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 5+sizes << " t \"max runtime\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 4+sizes << " t \"min runtime\" with linespoints\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"MFLOPS\"\nset output \"" << eps2name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 2+sizes << " t \"median FLOPS\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 8+sizes << " t \"mean FLOPS\" with linespoints\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"MB/s\"\nset output \"" << eps3name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 3+sizes << " t \"median transfer rate\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 2:" << 9+sizes << " t \"mean transfer rate\" with linespoints\n";
    }
    else if (plotcvx && not plotsvx)
    {
        ofs2 << "set terminal postscript eps color\nset key below\nset yrange [0:]\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"number of parts\"\nset ylabel \"time in sec.\"\nset output \"" << eps1name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 7+sizes << " t \"median runtime\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 5+sizes << " t \"max runtime\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 4+sizes << " t \"min runtime\" with linespoints\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"number of parts\"\nset ylabel \"MFLOPS\"\nset output \"" << eps2name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 2+sizes << " t \"median FLOPS\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 8+sizes << " t \"mean FLOPS\" with linespoints\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"number of parts\"\nset ylabel \"MB/s\"\nset output \"" << eps3name << "\"\n";
        ofs2 << "plot \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 3+sizes << " t \"median transfer rate\" with linespoints, \"BenchmarkPlotData\" index " << nr-1 << " using 1:" << 9+sizes << " t \"mean transfer rate\" with linespoints\n";
    }
    else if (plotcvx && plotsvx && (cnr == 1))
    {
        int xv(1), yv(1);
        std::list<std::string>::iterator ci = cores.begin();
        std::string st = *ci;
        ++ci;
        while (st != *ci)
        {
            ++xv;
            ++ci;
        }
        yv = int(cores.size() / xv);
        ofs2 << "set terminal postscript eps color\nset key below\nset contour base\nset surface\nset dgrid3d " << yv << "," << xv << ", \nshow contour\nset zrange [0:]\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"number of parts\"\nset ylabel \"Operand Size\"\nset zlabel \"time in sec.\"\nset output \"" << eps1name << "\"\n";
        ofs2 << "splot \"BenchmarkPlotData\" index " << nr-1 << " using 1:2:" << 7+sizes << " t \"median runtime\" with lines\n";
        ofs2 << "set title \"" << _id << "\"\nset zlabel \"MFLOPS\"\nset output \"" << eps2name << "\"\n";
        ofs2 << "splot \"BenchmarkPlotData\" index " << nr-1 << " using 1:2:" << 2+sizes << " t \"median FLOPS\" with lines\n";
        ofs2 << "set title \"" << _id << "\"\nset zlabel \"MB/s\"\nset output \"" << eps3name << "\"\n";
        ofs2 << "splot \"BenchmarkPlotData\" index " << nr-1 << " using 1:2:" << 3+sizes << " t \"median transfer rate\" with lines\n";
    }
    else if (plotcvx && plotsvx && (cnr > 1))
    {
        int noyellow(1);
        temp = *cores.begin();
        std::list<std::string>::iterator cit = cores.begin();
        ofs2 << "set terminal postscript eps color\nset key below\nset yrange [0:]\n";
        ofs2 << "set title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"time in sec.\"\nset output \"" << eps1name << "\"\n";
        for (int i(0) ; i < cnr ; ++i, ++noyellow)
        {
            if (noyellow == 6)
                ++noyellow;
            if (i != 0)
            {
                ofs2 << ", ";
                while (temp == *cit)
                    ++cit;
                temp = *cit;
            }
            else
                ofs2 << "plot ";
            ofs2 << "\"BenchmarkPlotData\" index " << nr - cnr + i << " using 2:" << 7+sizes << " t \"" + temp + " median runtime\" with linespoints " << noyellow;
            
        }
        temp = *cores.begin();
        cit = cores.begin();
        noyellow = 1;
        ofs2 << "\nset title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"MFLOPS\"\nset output \"" << eps2name << "\"\n";
        for (int i(0) ; i < cnr ; ++i, ++noyellow)
        {
            if (noyellow == 6)
                ++noyellow;
            if (i != 0)
            {
                ofs2 << ", ";
                while (temp == *cit)
                    ++cit;
                temp = *cit;
            }
            else
                ofs2 << "plot ";
            ofs2 << "\"BenchmarkPlotData\" index " << nr - cnr + i << " using 2:" << 2+sizes << " t \"" + temp + " median FLOPS\" with linespoints " << noyellow;
        }
        temp = *cores.begin();
        cit = cores.begin();
        noyellow = 1;
        ofs2 << "\nset title \"" << _id << "\"\nset xlabel \"Operand Size\"\nset ylabel \"MB/s\"\nset output \"" << eps3name << "\"\n";
        for (int i(0) ; i < cnr ; ++i, ++noyellow)
        {
            if (noyellow == 6)
                ++noyellow;
            if (i != 0)
            {
                ofs2 << ", ";
                while (temp == *cit)
                    ++cit;
                temp = *cit;
            }
            else
                ofs2 << "plot ";
            ofs2 << "\"BenchmarkPlotData\" index " << nr - cnr + i << " using 2:" << 3+sizes << " t \"" + temp + " median transfer rate\" with linespoints " << noyellow;
        }
    }
    else
    {
        std::cout << "Can't find a way to plot data automatically" << std::endl;
    }
    ofs2.close();
    ofs3 << "\t\t\\hline\n";
    ofs3 << "\t\\end{longtable}\n";
    ofs3 << "\t\\begin{figure}[H]\n";
    ofs3 << "\t\t\\includegraphics{" << eps1name << "}\n";
    ofs3 << "\t\\end{figure}\n";
    ofs3 << "\t\\begin{figure}[H]\n";
    ofs3 << "\t\t\\includegraphics{" << eps2name << "}\n";
    ofs3 << "\t\\end{figure}\n";
    ofs3 << "\t\\begin{figure}[H]\n";
    ofs3 << "\t\t\\includegraphics{" << eps3name << "}\n";
    ofs3 << "\t\\end{figure}\n";
    ofs3 << "\t\\end{center}\n";
    ofs3.close();
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

bool Benchmark::plots()
{
    return _plots;
}

int main(int argc, char** argv)
{
    int result=EXIT_SUCCESS;
    std::list<int> runrs;
    bool sse(true);
    bool cuda(true);
    bool cell(true);
    bool mc(true);
    bool sc(true);
    bool interface(false);
    bool plot(false);
    if ((argc == 2) && (honei::stringify(argv[1]) == "i"))
        interface = true;
    else if ((argc == 2) && (honei::stringify(argv[1]) == "plot"))
        plot = true;
    else if ((argc == 3) && ( ((honei::stringify(argv[1]) == "i") && (honei::stringify(argv[2]) == "plot")) || ((honei::stringify(argv[1]) == "plot") && (honei::stringify(argv[2]) == "i")) ) )
    {
        interface = true;
        plot = true;
    }
    else
    {
        if (argc > 1)
        {
            sse = false;
            cell = false;
            mc = false;
            sc = false;
            cuda = false;
            for(int i(1) ; i < argc ; ++i)
            {
                if (honei::stringify(argv[i]) == "sse")
                {
                    sse = true;
                }
                if (honei::stringify(argv[i]) == "cuda")
                {
                    cuda = true;
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
                if (honei::stringify(argv[i]) == "i")
                {
                    interface = true;
                }
                if (honei::stringify(argv[i]) == "plot")
                {
                    plot = true;
                }
            }
        }
    }
    for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_benchs()),i_end(BenchmarkList::instance()->end_benchs()) ; i != i_end ; )
    {
        if (sse && ((*i)->plots() == plot) && (((*i)->get_tag_name() == "sse") || ((*i)->get_tag_name() == "mc-sse")))
        {
            ++i;
            continue;
        }
        if (cell && ((*i)->plots() == plot) && ((*i)->get_tag_name() == "cell"))
        {
            ++i;
            continue;
        }
        if (cuda && ((*i)->plots() == plot) && ((*i)->get_tag_name() == "cuda"))
        {
            ++i;
            continue;
        }
            if (mc && ((*i)->plots() == plot) && ((*i)->get_tag_name() == "mc"))
        {
            ++i;
            continue;
        }
        if (sc && ((*i)->plots() == plot) && ((*i)->get_tag_name() == "cpu"))
        {
            ++i;
            continue;
        }
        i = BenchmarkList::instance()->erase(i);
    }
    if (BenchmarkList::instance()->begin_benchs() == BenchmarkList::instance()->end_benchs())
    {
        std::cout << "No relevant Benchmarks." << std::endl;
        return result;
    }
    if (interface)
    {
        int count = 1;
        std::cout << "Select Benchmark you'd like to add to runlist:" << std::endl;
        for (BenchmarkList::Iterator i(BenchmarkList::instance()->begin_benchs()),i_end(BenchmarkList::instance()->end_benchs()) ; i != i_end ; )
        {
            std::cout << count << ": " << (*i)->id() << std::endl;
            ++count;
            ++i;
        }
        std::cout << std::endl << "Choose: (1), (2), ..., (a)ll, (n)one" << std::endl;
        std::string a, tmp;
        int b;
        while (a != "n")
        {
            std::cout << "Selection: ";
            std::cin >> a;
            std::cout << "Added: ";
            if (a == "a")
            {
                std::cout << "all" << std::endl;
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
                        std::cout << b << " ";
                    }
                }
            }
            std::cout << std::endl;
            if (a != "n") std::cout << "Add more? ";
        }
    } 
    else 
    {
        for (int i = 0 ; i < BenchmarkList::instance()->bench_count() ; ++i)
        {
            runrs.push_back(i);
        }
    }
    if (plot)
    {
        time_t t;
        time(&t);
        std::ofstream ofs("RecentPlots.tex", std::ios_base::out);
        ofs << "\\documentclass{report}\n";
//        ofs << "\\usepackage{fullpage}\n";
        ofs << "\\usepackage{epsfig}\n";
        ofs << "\\usepackage{epstopdf}\n";
        ofs << "\\usepackage{float}\n";
        ofs << "\\usepackage{hyperref}\n";
        ofs << "\\usepackage{longtable}\n";
        ofs << "\\begin{document}\n";
        ofs << "\\tableofcontents\n";
        ofs << "\\newpage\n";
        ofs << "\t\\pagestyle{myheadings}\n"; 
        ofs << "\t\\markboth{}{" << ctime(&t) << "}\n";
        ofs.close();
    }
    int count = 0;
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
                    std::cout << (*i)->id() << ":" << std::endl;
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
    if (plot)
    {
        std::ofstream ofs1("RecentPlots.tex", std::ios_base::out | std::ios_base::app);
        ofs1 << "\\end{document}";
    }
    return result;
}
