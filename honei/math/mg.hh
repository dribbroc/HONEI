/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011-2012 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef MATH_GUARD_MG_HH
#define MATH_GUARD_MG_HH 1

#include <honei/math/operator.hh>
#include <honei/math/operator_list.hh>
#include <honei/util/profiler.hh>
#include <honei/math/methods.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/transposition.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/spai2.hh>

namespace honei
{
    struct MGDataIndex
    {
        public:
            static unsigned long internal_index_A(unsigned long base_index)
            {
                return base_index - 1;
            }
            static unsigned long base_index_A(unsigned long internal_index)
            {
                return internal_index + 1;
            }
            static unsigned long internal_index_Prol(unsigned long base_index)
            {
                return base_index - 2;
            }
            static unsigned long base_index_Prol(unsigned long internal_index)
            {
                return internal_index + 2;
            }
    };

    template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_, typename DataType_>
        struct MGData
        {
            public:
            std::vector<MatrixType_> A;
            std::vector<TransferContType_> resmat;
            std::vector<TransferContType_> prolmat;
            std::vector<PreconContType_> P;
            std::vector<VectorType_> b;
            std::vector<VectorType_> x;
            std::vector<VectorType_> d;
            std::vector<VectorType_> c;
            std::vector<VectorType_> store;
            std::vector<std::vector<VectorType_> > smoother_temp;
            unsigned long max_iters;
            unsigned long max_iters_coarse;
            unsigned long used_iters_coarse;
            unsigned long used_iters;
            unsigned long n_pre_smooth;
            unsigned long n_post_smooth;
            unsigned long min_level;
            DataType_ eps_relative;

            //storage used by cycle_control
            bool down_sweep;
            unsigned long max_iters_coarse_FMG;

            MGData(std::vector<MatrixType_> & systems,
                   std::vector<TransferContType_> & resmats,
                   std::vector<TransferContType_> & prolmats,
                   std::vector<PreconContType_> & precons,
                   std::vector<VectorType_> & rhss,
                   std::vector<VectorType_> & xs,
                   std::vector<VectorType_> & ds,
                   std::vector<VectorType_> & cs,
                   std::vector<VectorType_> & stores,
                   std::vector<std::vector<VectorType_> > & smoother_temps,
                   unsigned long p_max_iters,
                   unsigned long p_max_iters_coarse,
                   unsigned long p_n_pre_smooth,
                   unsigned long p_n_post_smooth,
                   unsigned long p_min_level,
                   DataType_ p_eps_rel
                   ) :
                A(systems),
                resmat(resmats),
                prolmat(prolmats),
                P(precons),
                b(rhss),
                x(xs),
                d(ds),
                c(cs),
                store(stores),
                smoother_temp(smoother_temps),
                max_iters(p_max_iters),
                max_iters_coarse(p_max_iters_coarse),
                used_iters_coarse(0),
                used_iters(0),
                n_pre_smooth(p_n_pre_smooth),
                n_post_smooth(p_n_post_smooth),
                min_level(p_min_level),
                eps_relative(p_eps_rel)
            {
                CONTEXT("When creating MGData:");
                ASSERT(p_min_level <= systems.size(), "Minimum level is larger then number of levels!");
            }


            template<typename MatrixSrc_, typename VectorSrc_, typename TransferSrc_, typename PreconSrc_, typename DataTypeSrc_>
            MGData(const MGData<MatrixSrc_, VectorSrc_, TransferSrc_, PreconSrc_, DataTypeSrc_> & other)
            {
                for (unsigned long i(0) ; i < other.A.size() ; ++i)
                {
                    MatrixType_ t1(other.A.at(i));
                    this->A.push_back(t1);
                }

                for (unsigned long i(0) ; i < other.resmat.size() ; ++i)
                {
                    TransferContType_ t2(other.resmat.at(i));
                    this->resmat.push_back(t2);
                }

                for (unsigned long i(0) ; i < other.prolmat.size() ; ++i)
                {
                    TransferContType_ t4(other.prolmat.at(i));
                    this->prolmat.push_back(t4);
                }

                for (unsigned long i(0) ; i < other.P.size() ; ++i)
                {
                    PreconContType_ t5(other.P.at(i));
                    this->P.push_back(t5);
                }

                for (unsigned long i(0) ; i < other.b.size() ; ++i)
                {
                    VectorType_ t1(other.b.at(i));
                    this->b.push_back(t1);
                    VectorType_ t2(other.x.at(i));
                    this->x.push_back(t2);
                    VectorType_ t3(other.d.at(i));
                    this->d.push_back(t3);
                    VectorType_ t4(other.c.at(i));
                    this->c.push_back(t4);
                    VectorType_ t5(other.store.at(i));
                    this->store.push_back(t5);
                }

                for(unsigned long i(0) ; i < other.smoother_temp.size() ; ++i)
                {
                    std::vector<VectorType_> temp;
                    for(unsigned long j(0) ; j < other.smoother_temp.at(i).size() ; ++j)
                    {
                        VectorType_ tv(other.smoother_temp.at(i).at(j));
                        temp.push_back(tv);
                    }
                    this->smoother_temp.push_back(temp);
                }

                this->max_iters = other.max_iters;
                this->max_iters_coarse = other.max_iters_coarse;
                this->n_pre_smooth = other.n_pre_smooth;
                this->n_post_smooth = other.n_post_smooth;
                this->min_level = other.min_level;
                this->eps_relative = other.eps_relative;
                this->used_iters = other.used_iters;
                this->used_iters_coarse = other.used_iters_coarse;
            }

            MGData copy()
            {
                std::vector<MatrixType_> Ar;
                std::vector<TransferContType_> resmatr;
                std::vector<TransferContType_> prolmatr;
                std::vector<PreconContType_> Pr;
                std::vector<VectorType_> br;
                std::vector<VectorType_> xr;
                std::vector<VectorType_> dr;
                std::vector<VectorType_> cr;
                std::vector<VectorType_> storer;
                std::vector<std::vector<VectorType_> > tempr;

                for (unsigned long i(0) ; i < this->A.size() ; ++i)
                {
                    Ar.push_back(this->A.at(i).copy());
                }

                for (unsigned long i(0) ; i < this->resmat.size() ; ++i)
                {
                    resmatr.push_back(this->resmat.at(i).copy());
                }

                for (unsigned long i(0) ; i < this->prolmat.size() ; ++i)
                {
                    prolmatr.push_back(this->prolmat.at(i).copy());
                }

                for (unsigned long i(0) ; i < this->P.size() ; ++i)
                {
                    Pr.push_back(this->P.at(i).copy());
                }

                for (unsigned long i(0) ; i < this->b.size() ; ++i)
                {
                    br.push_back(this->b.at(i).copy());
                    xr.push_back(this->x.at(i).copy());
                    dr.push_back(this->d.at(i).copy());
                    cr.push_back(this->c.at(i).copy());
                    storer.push_back(this->store.at(i).copy());
                }

                for(unsigned long i(0) ; i < this->smoother_temp.size() ; ++i)
                {
                    std::vector<VectorType_> temp(this->smoother_temp.at(i));
                    tempr.push_back(temp);
                }

                MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> result(Ar, resmatr, prolmatr, Pr, br, xr, cr, dr, storer, tempr, 1, 1000, 4, 4, 1, 1e-8);


                result.max_iters = this->max_iters;
                result.max_iters_coarse = this->max_iters_coarse;
                result.n_pre_smooth = this->n_pre_smooth;
                result.n_post_smooth = this->n_post_smooth;
                result.min_level = this->min_level;
                result.eps_relative = this->eps_relative;
                //result.used_iters = this->used_iters = other.used_iters;
                //this->used_iters_coarse = other.used_iters_coarse;

                return result;
            }

        };

    template<typename PreconContType_,
        typename MatrixType_,
        typename DataType_>
            struct PreconFill
            {
            };

    template<typename MatrixType_, typename DataType_>
        struct PreconFill<DenseVector<DataType_>, MatrixType_, DataType_>
        {
            public:

                static void value(unsigned long /*i*/,
                        std::vector<DenseVector<DataType_> > & target,
                        std::string /*filename*/,
                        std::string precon_suffix,
                        MatrixType_ & A,
                        DataType_ damping_factor)
                {
                    if (precon_suffix == "jac")
                    {
                        DenseVector<DataType_> current_P(A.rows(), damping_factor);
                        for(unsigned long j(0) ; j < current_P.size() ; ++j)
                        {
                            current_P[j] /= A(j, j);
                        }
                        target.push_back(current_P);
                    }
                    else if (precon_suffix == "spai")
                    {
                        DenseVector<DataType_> current_P(A.rows());
                        SparseMatrix<DataType_> As(A);
                        SPAI2<tags::CPU>::value_0(current_P, As);
                        Scale<tags::CPU>::value(current_P, damping_factor);
                        target.push_back(current_P);
                    }
                    else
                        throw InternalError("Preconditioner unknown!");
                }
        };

    template<typename DataType_>
        struct PreconFill<SparseMatrixELL<DataType_>, SparseMatrixELL<DataType_>, DataType_>
        {
            public:

                static void value(unsigned long i,
                        std::vector<SparseMatrixELL<DataType_> > & target,
                        std::string filename,
                        std::string precon_suffix,
                        SparseMatrixELL<DataType_> & /*A*/,
                        DataType_ damping_factor)
                {
                    std::string local_P_name(filename);
                    local_P_name += stringify(i+1);
                    local_P_name += "_";
                    local_P_name += precon_suffix;
                    local_P_name += ".ell";
                    SparseMatrixELL<DataType_> current_P(MatrixIO<io_formats::ELL>::read_matrix(local_P_name, DataType_(0)));
                    if (damping_factor != DataType_(1))
                    {
                        SparseMatrix<DataType_> sm_P(current_P);
                        Scale<tags::CPU>::value(sm_P, damping_factor);
                        SparseMatrixELL<DataType_> new_P(sm_P);
                        target.push_back(new_P);
                    }
                    else
                        target.push_back(current_P);
                }
        };

    template<typename Tag_, typename VectorType_>
    void remove_redundant_operators(OperatorList & ol)
    {
        for(unsigned long i(1) ; i < ol.size() ; ++i)
        {
            if(typeid(*ol[i]) == typeid(*ol[i - 1]))
                //if(!(typeid(*ol[i]) == typeid(FillOperator<Tag_, VectorType_>)))
                    ol.erase(i);
        }
    }

    void print_cycle(OperatorList & ol, unsigned long max_level, unsigned long min_level)
    {
        std::vector<unsigned long> transfers;
        for (unsigned long i(0) ; i < ol.size() ; ++i)
        {
            if (ol[i]->transfer_type() > 0)
                transfers.push_back(ol[i]->transfer_type());
        }

        DenseMatrix<unsigned long> graph((max_level - min_level + 1) * 2, transfers.size(), 0ul);
        unsigned long acx(0);
        for (unsigned long i(0) ; i < transfers.size() ; ++i)
        {
            switch (transfers.at(i))
            {
                case 1:
                    ++acx;
                    graph(acx, i) = transfers.at(i);
                    ++acx;
                    break;

                case 2:
                    --acx;
                    graph(acx, i) = transfers.at(i);
                    --acx;
                    break;

                case 3:
                    graph(acx, i) = transfers.at(i);
                    break;

                case 4:
                    graph(acx, i) = transfers.at(i);
                    break;
            }
        }

        std::cout<<"R: Restriction, P: Prolongation, S: Smoother, C: Solver"<<std::endl<<std::endl;
        for (unsigned long i(0) ; i < graph.rows() ; ++i)
        {
            for (unsigned long j(0) ; j < graph.columns() ; ++j)
            {
                switch(graph(i, j))
                {
                    case 0:
                        std::cout<<" ";
                        break;

                    case 1:
                        std::cout<<"R";
                        break;

                    case 2:
                        std::cout<<"P";
                        break;

                    case 3:
                        std::cout<<"S";
                        break;

                    case 4:
                        std::cout<<"C";
                        break;
                }
            }
            std::cout<<std::endl;
        }
    }

    template<typename Tag_,
        typename MatrixType_,
        typename VectorType_,
        typename TransferContType_,
        typename PreconContType_,
        typename MatIOType_,
        typename VecIOType_,
        typename DataType_>
            struct MGUtil
            {
                public:
                    static void configure(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & target, unsigned long max_iters,
                            unsigned long max_iters_coarse,
                            unsigned long n_pre_smooth,
                            unsigned long n_post_smooth,
                            unsigned long min_level,
                            DataType_ eps_relative)
                    {
                        target.max_iters = max_iters;
                        target.max_iters_coarse = max_iters_coarse;
                        target.n_pre_smooth = n_pre_smooth;
                        target.n_post_smooth = n_post_smooth;
                        target.min_level = min_level;
                        target.eps_relative = eps_relative;
                    }


                    static MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> load_data(std::string file_base, unsigned long max_level, DataType_ damping_factor, std::string precon_suffix)
                    {
                        CONTEXT("When creating MGData from file(s):");

                        std::vector<MatrixType_> A;
                        std::vector<TransferContType_> Prol;
                        std::vector<TransferContType_> Res;
                        std::vector<PreconContType_> P;
                        std::vector<VectorType_> b;
                        std::vector<VectorType_> x;
                        std::vector<VectorType_> d;
                        std::vector<VectorType_> c;
                        std::vector<VectorType_> store;
                        std::vector<std::vector<VectorType_> > smoother_temp;


                        std::string A_name(file_base);
                        std::string Prol_name(file_base);
                        std::string b_name(file_base);
                        std::string x_name(file_base);
                        std::string P_name(file_base);

                        A_name += "A_";

                        P_name += "A_";

                        Prol_name += "prol_";
                        b_name += "rhs_";
                        x_name += "init_";

                        for(unsigned long i(0) ; i <= max_level ; ++i)
                        {
                            if(i < max_level)
                            {
                                ///get system-matrix for level i
                                std::string local_A_name(A_name);
                                local_A_name += stringify(MGDataIndex::base_index_A(i));
                                local_A_name += ".ell";
                                MatrixType_ local_A(MatrixIO<MatIOType_>::read_matrix(local_A_name, DataType_(0)));
                                A.push_back(local_A);

                                PreconFill<PreconContType_, MatrixType_, DataType_>::value(i, P, P_name, precon_suffix, local_A, damping_factor);
                            }

                            ///get Prolmat P_{i}^{i+1}
                            if(i < max_level - 1)
                            {
                                std::string local_Prol_name(Prol_name);
                                local_Prol_name += stringify(MGDataIndex::base_index_Prol(i));
                                local_Prol_name += ".ell";
                                TransferContType_ local_Prol(MatrixIO<MatIOType_>::read_matrix(local_Prol_name, DataType_(0)));
                                Prol.push_back(local_Prol);

                                ///get Resmat R_{i+1}^{i} = (P_{i}^{i+1})^T
                                SparseMatrix<DataType_> local_preProl(Prol.at(i).copy());
                                SparseMatrix<DataType_> local_preRes(Prol.at(i).columns(), Prol.at(i).rows());
                                Transposition<Tag_>::value(local_preProl, local_preRes);
                                TransferContType_ local_Res(local_preRes);
                                Res.push_back(local_Res);
                            }

                            ///get vectors for level max_level
                            if(i == MGDataIndex::internal_index_A(max_level))
                            {
                                std::string local_b_name(b_name);
                                local_b_name += stringify(MGDataIndex::base_index_A(i));
                                VectorType_ max_b(VectorIO<VecIOType_>::read_vector(local_b_name, DataType_(0)));
                                b.push_back(max_b);

                                std::string local_x_name(x_name);
                                local_x_name += stringify(MGDataIndex::base_index_A(i));
                                VectorType_ max_x(VectorIO<VecIOType_>::read_vector(local_x_name, DataType_(0)));
                                //VectorType_ max_x(max_b.size(), DataType_(0));
                                x.push_back(max_x);

                                VectorType_ zero(A.at(i).rows(), DataType_(0));
                                d.push_back(zero.copy());
                                ///Set c = x on finest level
                                //c.push_back(max_x.copy());
                                c.push_back(zero.copy());
                                store.push_back(zero.copy());
                                std::vector<VectorType_> smtl;
                                smoother_temp.push_back(smtl);
                            }
                            else if(i == MGDataIndex::internal_index_A(1))
                            {
                                std::string local_b_name(b_name);
                                local_b_name += stringify(MGDataIndex::base_index_A(i));
                                VectorType_ max_b(VectorIO<VecIOType_>::read_vector(local_b_name, DataType_(0)));
                                b.push_back(max_b);

                                std::string local_x_name(x_name);
                                local_x_name += stringify(MGDataIndex::base_index_A(i));
                                VectorType_ max_x(VectorIO<VecIOType_>::read_vector(local_x_name, DataType_(0)));
                                //VectorType_ max_x(max_b.size(), DataType_(0));
                                x.push_back(max_x);

                                VectorType_ zero(A.at(i).rows(), DataType_(0));
                                d.push_back(zero.copy());
                                ///Set c = x on finest level
                                //c.push_back(max_x.copy());
                                c.push_back(zero.copy());
                                store.push_back(zero.copy());
                                std::vector<VectorType_> smtl;
                                smoother_temp.push_back(smtl);
                            }
                            else if(i < MGDataIndex::internal_index_A(max_level))
                            {
                                ///get vectors for level i
                                VectorType_ zero(A.at(i).rows(), DataType_(0));
                                b.push_back(zero.copy());
                                x.push_back(zero.copy());
                                d.push_back(zero.copy());
                                c.push_back(zero.copy());
                                store.push_back(zero.copy());
                                std::vector<VectorType_> smtl;
                                smoother_temp.push_back(smtl);
                            }

                            if(typeid(PreconContType_) == typeid(MatrixType_))
                            {
                                std::string P_name(file_base);
                                P_name += "A_spai_"; //TODO: we should rename all precalculated precons to "P_"; what if there are more types?
                                //TODO
                            }
                        }

                        MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> result(A, Res, Prol, P, b, x, d, c, store, smoother_temp, 0, 0, 0, 0, 0, DataType_(0.));
                        return result;
                    }
            };

    //TODO: think of more flexible way (more than one smoothertype, ...) ; later, assembly is done here!
    template<typename Tag_,
        typename CycleShape_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation
            {
            };

    //specialise by cycle-shape
    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::V::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }


            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        _build_cycle(data.x, data.b, data.A.size(), cycle, data);
                        return cycle;
                    }
        };

    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::W::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            if (level - 1 != data.min_level)
                                _build_cycle(data.c, data.d, level - 1, cycle, data);
                            _build_cycle(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        _build_cycle(data.x, data.b, data.A.size(), cycle, data);
                        //remove_redundant_operators<Tag_, VectorType_>(cycle);
                        return cycle;
                    }
        };

    //specialise by cycle-shape
    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::F::V::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_V(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_V(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_F(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_F(data.c, data.d, level - 1, cycle, data);
                            _build_cycle_V(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        _build_cycle_F(data.x, data.b, data.A.size(), cycle, data);
                        //remove_redundant_operators<Tag_, VectorType_>(cycle);
                        return cycle;
                    }
        };

    //specialise by cycle-shape
    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::F::W::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_W(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            ////std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            if (level - 1 != data.min_level)
                                _build_cycle_W(data.c, data.d, level - 1, cycle, data);
                            _build_cycle_W(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol((level)))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            ////std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_F(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.smoother_temp.at(MGDataIndex::internal_index_A(level)).at(0), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_F(data.c, data.d, level - 1, cycle, data);
                            _build_cycle_W(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        _build_cycle_F(data.x, data.b, data.A.size(), cycle, data);
                        //remove_redundant_operators<Tag_, VectorType_>(cycle);
                        return cycle;
                    }
        };

    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::FMG::V::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_V(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_V(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_F(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_F(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                        }
                        _build_cycle_V(x, b, level, cycle, data);
                    }
            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        data.down_sweep = true;
                        data.max_iters_coarse_FMG = 1000;
                        _build_cycle_F(data.x, data.b, data.A.size(), cycle, data);
                        //remove_redundant_operators<Tag_, VectorType_>(cycle);
                        return cycle;
                    }
        };

    template<typename Tag_,
        typename CoarseGridSolverType_,
        typename SmootherType_,
        typename ResType_,
        typename ProlType_,
        typename DataType_>
            struct MGCycleCreation<Tag_,
        methods::CYCLE::FMG::W::STATIC,
        CoarseGridSolverType_,
        SmootherType_,
        ResType_,
        ProlType_,
        DataType_>
        {
            private:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_W(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Presmoothing
                            //std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_pre_smooth) );


                            ///Defect
                            ////std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            if (level - 1 != data.min_level)
                                _build_cycle_W(data.c, data.d, level - 1, cycle, data);
                            _build_cycle_W(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol((level)))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                            ///Postsmoothing
                            ////std::cout << "Smoother Accessing " << level << std::endl;
                            cycle.push_back(new SmootherOperator<SmootherType_, MatrixType_, VectorType_, PreconContType_>(
                                        data.A.at(MGDataIndex::internal_index_A(level)),
                                        data.P.at(MGDataIndex::internal_index_A(level)),
                                        b.at(MGDataIndex::internal_index_A(level)),
                                        x.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_0.at(MGDataIndex::internal_index_A(level)),
                                        //data.temp_1.at(MGDataIndex::internal_index_A(level)),
                                        data.smoother_temp.at(MGDataIndex::internal_index_A(level)),
                                        data.n_post_smooth) );
                        }
                    }

                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static void _build_cycle_F(
                            std::vector<VectorType_> & x,
                            std::vector<VectorType_> & b,
                            unsigned long level,
                            OperatorList & cycle,
                            MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        if(level == (data.min_level))
                        {
                            //std::cout << "Solver Accessing " << data.min_level << std::endl;
                            cycle.push_back(new SolverOperator<CoarseGridSolverType_, MatrixType_, VectorType_, PreconContType_, DataType_>(
                                        data.A.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.P.at(MGDataIndex::internal_index_A(data.min_level)),
                                        b.at(MGDataIndex::internal_index_A(data.min_level)),
                                        x.at(MGDataIndex::internal_index_A(data.min_level)),
                                        data.max_iters_coarse,
                                        data.used_iters_coarse,
                                        data.eps_relative) );
                        }
                        else
                        {
                            ///Defect
                            //std::cout << "Defect Accessing " << level << std::endl;
                            cycle.push_back(new DefectOperator<Tag_, MatrixType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), b.at(MGDataIndex::internal_index_A(level)), data.A.at(MGDataIndex::internal_index_A(level)), x.at(MGDataIndex::internal_index_A(level))));
                            ///Restriction
                            //std::cout << " Restrict Accessing " << level << std::endl;
                            //std::cout << " Restrict Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ResType_, TransferContType_, VectorType_>(data.d.at(MGDataIndex::internal_index_A(level - 1)) , data.store.at(MGDataIndex::internal_index_A(level)), data.resmat.at(MGDataIndex::internal_index_Prol(level))));

                            ///Recursion
                            ///all vectors in c have to be initialised with 0
                            cycle.push_back(new FillOperator<Tag_, VectorType_>(data.c.at(MGDataIndex::internal_index_A(level - 1))));
                            _build_cycle_F(data.c, data.d, level - 1, cycle, data);

                            ///Prolongation
                            //std::cout << "Prol Accessing " << level << std::endl;
                            //std::cout << "Prol Accessing " << level - 1 << std::endl;
                            cycle.push_back(new TransferOperator<ProlType_, TransferContType_, VectorType_>(data.store.at(MGDataIndex::internal_index_A(level)), data.c.at(MGDataIndex::internal_index_A(level - 1)), data.prolmat.at(MGDataIndex::internal_index_Prol(level))));
                            //std::cout << "Sum Accessing " << level << std::endl;
                            cycle.push_back(new SumOperator<Tag_, VectorType_>(x.at(MGDataIndex::internal_index_A(level)), data.store.at(MGDataIndex::internal_index_A(level))));

                        }
                        _build_cycle_W(x, b, level, cycle, data);
                    }

            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_>
                    static OperatorList value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data)
                    {
                        CONTEXT("When evaluating MGCycleCreation:");

                        OperatorList cycle;
                        data.down_sweep = true;
                        data.max_iters_coarse_FMG = 1000;
                        _build_cycle_F(data.x, data.b, data.A.size(), cycle, data);
                        //remove_redundant_operators<Tag_, VectorType_>(cycle);
                        return cycle;
                    }
        };



    template<typename Tag_, typename NormType_>
        struct MGSolver
        {
            public:
                template<typename MatrixType_, typename VectorType_, typename TransferContType_, typename PreconContType_, typename DataType_>
                    static void value(MGData<MatrixType_, VectorType_, TransferContType_, PreconContType_, DataType_> & data, OperatorList & cycle)
                    {
                        CONTEXT("When solving linear system with MG :");
                        ASSERT(cycle.size() > 0, "OperatorList is empty!");
                        PROFILER_START("MGSolver");

                        VectorType_ r(data.x.at(data.x.size() - 1).size());
                        Defect<Tag_>::value(r, data.b.at(data.b.size() - 1), data.A.at(data.A.size() - 1), data.x.at(data.x.size() - 1));
                        DataType_ rnorm_initial = NormType_::value(r);
                        std::cout << "Initial norm: " << rnorm_initial << std::endl;
                        DataType_ rnorm_current(1e16);

                        //std::cout << "starting cycles" << std::endl;
                        for(unsigned long i(0) ; i < data.max_iters ; ++i)
                        {
                            cycle.value();
                            Defect<Tag_>::value(r, data.b.at(data.b.size() - 1), data.A.at(data.A.size() - 1), data.x.at(data.x.size() - 1));
                            rnorm_current = NormType_::value(r);
                            //std::cout << "DEFECTNORM: " << rnorm_current << std::endl;

                            data.used_iters = i + 1;

                            if(rnorm_current < data.eps_relative * rnorm_initial)
                                break;
                        }

                        PROFILER_STOP("MGSolver");
                        std::cout << "Final norm: " << rnorm_current << std::endl;
                    }
        };

    struct MGSmoother
    {
        public:
            static void value(OperatorList & cycle, unsigned long max_iters)
            {
                CONTEXT("When smoothing linear system with MG :");
                ASSERT(cycle.size() > 0, "OperatorList is empty!");
                PROFILER_START("MGSmoother");

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cycle.value();
                }
                PROFILER_STOP("MGSmoother");
            }
    };
}

#endif
