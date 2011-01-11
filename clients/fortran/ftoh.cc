/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <iostream>
#include <honei/la/dense_vector.hh>
#include <honei/la/sum.hh>
#include <honei/math/multigrid.hh>

using namespace honei;

extern "C"
{
    void sum_(double * a, double * b, long * size);
    void is_smoother_(bool* flag);
    void set_border_mask_(int * mask);
    void set_lvl_range_(int* min, int* max);
    void set_n_max_iter_(int* iter);
    void set_initial_zero_(bool* initz);
    void set_tolerance_(double* eps);
    void set_convergence_check_(bool* convc);
    void set_n_pre_smooth_(int* nps);
    void set_n_post_smooth_(int* nps);
    void set_n_max_iter_coarse_(int* iter);
    void set_tolerance_coarse_(double* eps);
    void set_adapt_correction_factor_(double *factor);
    void push_c_(double* data, long* size);
    void push_d_(double* data, long* size);
    void push_rhs_(double* data, long* size);
    void push_x_(double* data, long* size);
    void push_diags_inverted_(double* data, long* size);
    void push_a_(long* rows, long* columns, long* nz,
            long* row_indices, long* column_indices, double* data);
    void push_prolmats_(long* rows, long* columns, long* nz,
            long* row_indices, long* column_indices, double* data);
    void push_resmats_(long* rows, long* columns, long* nz,
            long* row_indices, long* column_indices, double* data);


    void run_();
}

void sum_(double * a, double * b, long * size)
{
    //DenseVector<double> av(*size);
    //DenseVector<double> bv(*size);
    DenseVector<float> av(*size);
    DenseVector<float> bv(*size);
    for (int i(0) ; i < *size ; ++i)
    {
        av[i] = a[i];
        bv[i] = b[i];
    }
    Sum<tags::GPU::CUDA>::value(av, bv);
    av.lock(lm_read_only);
    av.unlock(lm_read_only);
    for (int i(0) ; i < *size ; ++i)
    {
        a[i] = av[i];
    }
}

class MG
{
    private:

        MG()
        {
        }

    public:
        /*static MG* instance()
        {
            static MG result;
            return &result;
        }*/
        static MGInfo<double, SparseMatrixELL<double> >* info()
        {
            static MGInfo<double, SparseMatrixELL<double> > info;
            return &info;
        }
};

void is_smoother_(bool* flag)
{
    MG::info()->is_smoother = *flag;
}

void set_border_mask_(int* mask)
{
    for (unsigned long i(0) ; i < 8 ; ++i)
    {
        (MG::info()->macro_border_mask)[i] = mask[i];
    }
}

void set_lvl_range_(int* min, int* max)
{
    MG::info()->min_level = *min;
    MG::info()->max_level = *max;
}

void set_n_max_iter_(int* iter)
{
    MG::info()->n_max_iter = *iter;
}

void set_initial_zero_(bool* initz)
{
    MG::info()->initial_zero = *initz;
}

void set_tolerance_(double* eps)
{
    MG::info()->tolerance = *eps;
}

void set_convergence_check_(bool* convc)
{
    MG::info()->convergence_check = *convc;
}

void set_n_pre_smooth_(int* nps)
{
    MG::info()->n_pre_smooth = *nps;
}

void set_n_post_smooth_(int* nps)
{
    MG::info()->n_pre_smooth = *nps;
}

void set_n_max_iter_coarse_(int* iter)
{
    MG::info()->n_max_iter_coarse = *iter;
}

void set_tolerance_coarse_(double* eps)
{
    MG::info()->tolerance_coarse = *eps;
}

void set_adapt_correction_factor_(double *factor)
{
    MG::info()->adapt_correction_factor = *factor;
}

void push_c_(double* data, long* size)
{
    DenseVector<double> c(*size);
    for (long i(0) ; i < *size ; ++i)
    {
        c[i] = data[i];
    }
    MG::info()->c.push_back(c);
}

void push_d_(double* data, long* size)
{
    DenseVector<double> d(*size);
    for (long i(0) ; i < *size ; ++i)
    {
        d[i] = data[i];
    }
    MG::info()->d.push_back(d);
}

void push_rhs_(double* data, long* size)
{
    DenseVector<double> rhs(*size);
    for (long i(0) ; i < *size ; ++i)
    {
        rhs[i] = data[i];
    }
    MG::info()->rhs.push_back(rhs);
}

void push_x_(double* data, long* size)
{
    DenseVector<double> x(*size);
    for (long i(0) ; i < *size ; ++i)
    {
        x[i] = data[i];
    }
    MG::info()->x.push_back(x);
}

void push_diags_inverted_(double* data, long* size)
{
    DenseVector<double> diags_inverted(*size);
    for (long i(0) ; i < *size ; ++i)
    {
        diags_inverted[i] = data[i];
    }
    MG::info()->diags_inverted.push_back(diags_inverted);
}

void push_a_(long* rows, long* columns, long* nz,
        long* row_indices, long* column_indices, double* data)
{
    DenseVector<unsigned long> row_dv(*nz);
    DenseVector<unsigned long> column_dv(*nz);
    DenseVector<double> data_dv(*nz);
    for (long i(0) ; i < *nz ; ++i)
    {
        row_dv[i] = row_indices[i];
        column_dv[i] = column_indices[i];
        data_dv[i] = data[i];
    }

    SparseMatrix<double> sm(*rows, *columns, row_dv.elements(), column_dv.elements(), data_dv.elements(), *nz);
    SparseMatrixELL<double> smell(sm);
    MG::info()->a.push_back(smell);
}

void push_prolmats_(long* rows, long* columns, long* nz,
        long* row_indices, long* column_indices, double* data)
{
    DenseVector<unsigned long> row_dv(*nz);
    DenseVector<unsigned long> column_dv(*nz);
    DenseVector<double> data_dv(*nz);
    for (long i(0) ; i < *nz ; ++i)
    {
        row_dv[i] = row_indices[i];
        column_dv[i] = column_indices[i];
        data_dv[i] = data[i];
    }

    SparseMatrix<double> sm(*rows, *columns, row_dv.elements(), column_dv.elements(), data_dv.elements(), *nz);
    SparseMatrixELL<double> smell(sm);
    MG::info()->prolmats.push_back(smell);
}

void push_resmats_(long* rows, long* columns, long* nz,
        long* row_indices, long* column_indices, double* data)
{
    DenseVector<unsigned long> row_dv(*nz);
    DenseVector<unsigned long> column_dv(*nz);
    DenseVector<double> data_dv(*nz);
    for (long i(0) ; i < *nz ; ++i)
    {
        row_dv[i] = row_indices[i];
        column_dv[i] = column_indices[i];
        data_dv[i] = data[i];
    }

    SparseMatrix<double> sm(*rows, *columns, row_dv.elements(), column_dv.elements(), data_dv.elements(), *nz);
    SparseMatrixELL<double> smell(sm);
    MG::info()->resmats.push_back(smell);
}




void run_()
{
    DenseVector<double> rhs(MG::info()->rhs[MG::info()->max_level]);
    DenseVector<double> result(rhs.size(), double(0));
    SparseMatrixELL<double> system(MG::info()->a[MG::info()->max_level]);
    Multigrid<tags::CPU::SSE, tags::CPU::SSE, methods::NONE, methods::JAC, methods::CYCLE::V, methods::FIXED >::value(system, rhs, result, (unsigned long)11, std::numeric_limits<double>::epsilon(), *MG::info());
}
