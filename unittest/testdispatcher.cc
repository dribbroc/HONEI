/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "testdispatcher.hh"
#include "tests.hh"
using namespace std;

bool TestDispatcher::start_test(int typ, BaseTest * instances[], int instances_size)
{
    bool result = false;
    switch (typ)
    {
        case Tests::all: //instructions
                        break;
                        
        case Tests::vector_rw: /*instructions
                        Vector vector = new Vector (0);
                        _vector_rw(instances, instances_size, vector, 100)
                        */
                        break;
                        
        case Tests::matrix_rw: //instructions
                        break;                                            

        default:        //instructions
                        break;
    }

    return result;
}

/*bool TestDispatcher::_vector_rw(BaseTest * instances[],int instances_size, Vector vector , int operation_count);
{
    // operationen erzeugen und ausführen,ergebnis muss ursprungsvector sein
    unsigned long vsize= vector.size();
    Vector vresult;
    double value;
    int position;
    srand((unsigned)time(0));
    
    for (int j=0 , j<instance_size , ++j)
    {   
        vresult=vector;
        for (int i=0 , i<operation_count , ++i)
        {   
            operator=rand();
            position=rand()%vsize;
            vresult = instances[j].vector_rw(vector, position, value);
            if (vresult[position] != value) return false;
        }
    }
    return true;
}*/

/*Vector * TestDispatcher::_random_dense_vector(unsigned long size);
{
    Vector result = new Vector;
    srand((unsigned)time(0));
    for (unsigned long i=0 ; i<size ; ++i)
    {
        result[i]=rand();
    }
}*/
/*Vector * TestDispatcher::_random_dense_matrix(unsigned long columns, unsigned long rows);
{
    Matrix result = new Matrix;
    srand((unsigned)time(0));
    for (unsigned long i=0 ; i<columns ; ++i)
    {
        for (unsigned long j=0 ; j<rows ; ++j)
        {
            result[i,j]=rand();
        }
    }
}*/
