/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "testdispatcher.hh"
#include "tests.hh"
using namespace std;

bool TestDispatcher::start_test(int typ, BaseTest * instances[], int instances_size)
{
    switch (typ)
    {
        case Tests::all: //instructions
                        break;

        default:        //instructions
                        break;
    }

    return false;
}
/*Vector * TestDispatcher::_random_dense_vector(unsigned long size);
{
    Vector result = new Vector;
    srand((unsigned)time(0));
    for (unsigned long i=0 ; i<size ; ++i)
    {
        result[i]=rand();
    }
}
*/
