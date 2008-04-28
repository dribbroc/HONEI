/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBSWE_GUARD_ASSEMBLY_PROCESSING_HH
#define LIBSWE_GUARD_ASSEMBLY_PROCESSING_HH 1

#include <honei/util/attributes.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/product.hh>
#include <honei/swe/directions.hh>
#include <honei/swe/limiter.hh>
/**
 * \file
 * Implementation of assembly processing functions for RelaxSolver.
 *
 * \ingroup grplibswe
 **/

namespace honei
{
    namespace assembly_types
    {
        class MAIN
        {
            class M1M3;
            class M2M4;
        };
        class QUICK
        {
            class M6;
            class M8;
        };
    }

    template <typename Tag_, typename Type_>
    struct AssemblyProcessing
    {
    };

    /**
     * \brief Implementation of source processing.
     *
     * \ingroup grplibswe
     *
     **/
    template <typename Tag_>
    struct AssemblyProcessing<Tag_, assembly_types::MAIN::M1M3>
    {

        public:
            template <typename WorkPrec_>
            static inline void value(BandedMatrix<WorkPrec_> &  m1, BandedMatrix<WorkPrec_> & m3, DenseVector<WorkPrec_> & u, DenseVector<WorkPrec_> & v, WorkPrec_ delta_t, WorkPrec_ delta_x, unsigned long d_width, unsigned long d_height, DenseVector<WorkPrec_> & c)
            {
                CONTEXT("When processing RelaxSolver assembly, main stage, matrices m_1, m_3:");

                ///The bands containing data.
                DenseVector<WorkPrec_> m1diag(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m1bandPlus1(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m1bandPlus2(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m1bandMinus1(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m3diag(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m3bandPlus1(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m3bandPlus2(u.size(), WorkPrec_(0));
                DenseVector<WorkPrec_> m3bandMinus1(u.size(), WorkPrec_(0));

                WorkPrec_ * m1diag_e = m1diag.elements();
                WorkPrec_ * m1bandPlus1_e = m1bandPlus1.elements();
                WorkPrec_ * m1bandPlus2_e = m1bandPlus2.elements();
                WorkPrec_ * m1bandMinus1_e = m1bandMinus1.elements();
                WorkPrec_ * m3diag_e = m3diag.elements();
                WorkPrec_ * m3bandPlus1_e = m3bandPlus1.elements();
                WorkPrec_ * m3bandPlus2_e = m3bandPlus2.elements();
                WorkPrec_ * m3bandMinus1_e = m3bandMinus1.elements();

                WorkPrec_ * u_e = u.elements();
                WorkPrec_ * v_e = v.elements();
                WorkPrec_ * c_e = c.elements();

                ///Necessary values to be temporarily saved.
                WorkPrec_ HONEI_ALIGNED(16) tempPlus[(unsigned long)(3)];
                WorkPrec_ HONEI_ALIGNED(16) tempTopPlus[(unsigned long)(3)];

                WorkPrec_ HONEI_ALIGNED(16) phiPlusOld[(unsigned long)(3)];
                WorkPrec_ HONEI_ALIGNED(16) phiPlusNew[(unsigned long)(3)];
                WorkPrec_ HONEI_ALIGNED(16) phiMinusNew[(unsigned long)(3)];

                WorkPrec_ HONEI_ALIGNED(16) tempMinus[(unsigned long)(3)];
                WorkPrec_ HONEI_ALIGNED(16) tempTopMinus[(unsigned long)(3)];

                for(int i(0); i < 3; ++i)
                {
                    tempPlus[i] = WorkPrec_(0);
                    tempTopPlus[i] = WorkPrec_(0);
                    phiPlusOld[i] = WorkPrec_(0);
                    phiPlusNew[i] = WorkPrec_(0);
                    phiMinusNew[i] = WorkPrec_(0);
                    tempMinus[i] = WorkPrec_(0);
                    tempTopMinus[i] = WorkPrec_(0);
                }

                WorkPrec_ HONEI_ALIGNED(16)phiMinusOld;
                WorkPrec_ HONEI_ALIGNED(16)temp;
                WorkPrec_ HONEI_ALIGNED(16)tempTop;
                WorkPrec_ HONEI_ALIGNED(16)prefac = delta_t/(4*delta_x);

                ///Needed Iterators.
                //typename DenseVector<WorkPrec_>::ElementIterator d(m1diag.begin_elements());
                unsigned long d(0);
                //typename DenseVector<WorkPrec_>::ConstElementIterator i(u.begin_elements());
                unsigned long i(0);

                //typename DenseVector<WorkPrec_>::ElementIterator b1(m1bandPlus1.begin_elements());
                unsigned long b1(0);

                //typename DenseVector<WorkPrec_>::ElementIterator b2(m1bandPlus2.begin_elements());
                unsigned long b2(0);

                //typename DenseVector<WorkPrec_>::ElementIterator bminus1(m1bandMinus1.begin_elements());

                unsigned long bminus1(0);

                /*
                ///Iterate through the vectors in order to avoid boundary access.
                for( ; i.index() < 6*(d_width+4); ++i);
                for( ; d.index() < 6*(d_width+4); ++d);
                for( ; b1.index() < 6*(d_width+4); ++b1);
                for( ; b2.index() < 6*(d_width+4); ++b2);
                for( ; bminus1.index() < 6*(d_width+4); ++bminus1);
                */
                i += 6*(d_width+4);
                d += 6*(d_width+4);
                b1 += 6*(d_width+4);
                b2 += 6*(d_width+4);
                bminus1 += 6*(d_width+4);

                while(i < 3*(d_width + 4) * (d_height + 2))
                {
                    for(unsigned long k=0; k<3; ++k)
                    {
                        tempPlus[k]= (v_e)[i] + (c_e)[k]*(u_e)[i];
                        ++i;++d;++b1;++b2;++bminus1;
                    }

                    for(unsigned long k=0; k<3; ++k)
                    {
                        temp= (v_e)[i] + (c_e)[k]*(u_e)[i];
                        tempTopPlus[k] = temp - tempPlus[k];
                        tempPlus[k] = temp;
                        tempMinus[k] =  (v_e)[i] - (c_e)[k]*(u_e)[i];
                        ++i;++d;++b1;++b2;++bminus1;
                    }

                    for(unsigned long k=0; k<3; ++k)
                    {
                        temp= (v_e)[i] - (c_e)[k]*(u_e)[i];
                        tempTopMinus[k] = temp - tempMinus[k];
                        tempMinus[k] = temp;
                        temp= (v_e)[i] + (c_e)[k]*(u_e)[i];
                        tempTop = temp - tempPlus[k];
                        if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon() )
                        {
                            phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                        }
                        else
                        {
                            phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                        }
                        tempPlus[k]=temp;
                        tempTopPlus[k]=tempTop;
                        ++i;
                    }

                    for(unsigned long k=0; k<3; ++k)
                    {
                        temp = (v_e)[i] - (c_e)[k]*((u_e)[i]);
                        tempTop = temp - tempMinus[k];
                        if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon() )
                        {
                            phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);
                        }
                        else
                        {
                            phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/std::numeric_limits<WorkPrec_>::epsilon());
                        }
                        tempMinus[k]=temp;
                        tempTopMinus[k] = tempTop;
                        temp  = (v_e)[i] + (c_e)[k]* (u_e)[i];
                        tempTop = temp - tempPlus[k];
                        if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon() )
                        {
                            phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                        }
                        else
                        {
                            phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                        }
                        tempPlus[k]= temp;
                        tempTopPlus[k] = tempTop;
                        ++i;
                    }

                    for(unsigned long x = 0; x < d_width; ++x)
                    {
                        for(unsigned long k =0; k<3; ++k)
                        {
                            temp = prefac *(2 - phiPlusOld[k]);
                            m1bandMinus1_e[bminus1] =temp;
                            m3bandMinus1_e[bminus1] =temp * (c_e)[k];
                            m1diag_e[d] = prefac * (phiPlusNew[k] + phiPlusOld[k] + phiMinusNew[k]);
                            m3diag_e[d] = -(c_e)[k] * prefac *(4 - phiPlusNew[k] - phiPlusOld[k] + phiMinusNew[k]);

                            phiPlusOld[k]= phiPlusNew[k];
                            phiMinusOld = phiMinusNew[k];
                            temp = (v_e)[i] - (c)[k]*(u_e)[i];
                            tempTop = temp - tempMinus[k];

                            if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon() )
                            {
                                phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);
                            }
                            else
                            {
                                phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/std::numeric_limits<WorkPrec_>::epsilon());
                            }

                            tempMinus[k]=temp;
                            tempTopMinus[k] = tempTop;
                            temp  = (v_e)[i] + (c_e)[k]* (u_e)[i];
                            tempTop = temp - tempPlus[k];
                            if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon() )
                            {
                                phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                            }
                            else
                            {
                                phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());
                            }

                            tempPlus[k]= temp;
                            tempTopPlus[k] = tempTop;
                            m1bandPlus1_e[b1] = prefac * (-2 - phiPlusOld[k] - phiMinusOld - phiMinusNew[k]);
                            m3bandPlus1_e[b1] = (c_e)[k]*prefac * (2 - phiPlusOld[k] + phiMinusOld + phiMinusNew[k]);
                            m1bandPlus2_e[b2] = prefac* phiMinusNew[k];
                            m3bandPlus2_e[b2] = (c_e)[k] * prefac *(-phiMinusNew[k]);
                            ++i;++d;++b1;++b2;++bminus1;
                        }
                    }
                    // Skip ghost-cells
                    ++d;++b1;++b2;++bminus1;
                    ++d;++b1;++b2;++bminus1;
                    ++d;++b1;++b2;++bminus1;
                    ++d;++b1;++b2;++bminus1;
                    ++d;++b1;++b2;++bminus1;
                    ++d;++b1;++b2;++bminus1;
                }

                m1.insert_band(0, m1diag);
                m1.insert_band(3, m1bandPlus1);
                m1.insert_band(6, m1bandPlus2);
                m1.insert_band((-3), m1bandMinus1);
                m3.insert_band(0, m3diag);
                m3.insert_band(3, m3bandPlus1);
                m3.insert_band(6, m3bandPlus2);
                m3.insert_band((-3), m3bandMinus1);
#ifdef SOLVER_VERBOSE
                cout << "M_1:" << stringify(m1.band(ulint(0))) << endl;
                cout << "M_1:" << stringify(m1.band(ulint(3))) << endl;
                cout << "M_1:" << stringify(m1.band(ulint(6))) << endl;
                cout << "M_1:" << stringify(m1.band(ulint(-3))) << endl;
                std::cout << "Finished Matrix Assembly 1.\n";
#endif

            }
    };

    template <typename Tag_>
    struct AssemblyProcessing<Tag_, assembly_types::MAIN::M2M4>
    {

        public:
            template <typename WorkPrec_>
            static inline void value(BandedMatrix<WorkPrec_> &  m2, BandedMatrix<WorkPrec_> & m4, DenseVector<WorkPrec_> & u, DenseVector<WorkPrec_> & w, WorkPrec_ delta_t, WorkPrec_ delta_y, unsigned long d_width, unsigned long d_height, DenseVector<WorkPrec_> & dv)
            {
            CONTEXT("When processing RelaxSolver assembly, main stage, matrices m_2, m_4:");

            ///The bands containing data.
            DenseVector<WorkPrec_> m2diag(u.size(), WorkPrec_(0));      //zero
            DenseVector<WorkPrec_> m2bandPlus1(u.size(), WorkPrec_(0)); //one
            DenseVector<WorkPrec_> m2bandPlus2(u.size(), WorkPrec_(0)); //two
            DenseVector<WorkPrec_> m2bandMinus1(u.size(), WorkPrec_(0));//three
            DenseVector<WorkPrec_> m4diag(u.size(), WorkPrec_(0));      //zero
            DenseVector<WorkPrec_> m4bandPlus1(u.size(), WorkPrec_(0)); //one
            DenseVector<WorkPrec_> m4bandPlus2(u.size(), WorkPrec_(0)); //two
            DenseVector<WorkPrec_> m4bandMinus1(u.size(), WorkPrec_(0));//three

            WorkPrec_ * m2diag_e = m2diag.elements();
            WorkPrec_ * m2bandPlus1_e = m2bandPlus1.elements();
            WorkPrec_ * m2bandPlus2_e = m2bandPlus2.elements();
            WorkPrec_ * m2bandMinus1_e = m2bandMinus1.elements();
            WorkPrec_ * m4diag_e = m4diag.elements();
            WorkPrec_ * m4bandPlus1_e = m4bandPlus1.elements();
            WorkPrec_ * m4bandPlus2_e = m4bandPlus2.elements();
            WorkPrec_ * m4bandMinus1_e = m4bandMinus1.elements();

            WorkPrec_ * u_e = u.elements();
            WorkPrec_ * w_e = w.elements();
            WorkPrec_ * d_e = dv.elements();

            WorkPrec_ HONEI_ALIGNED(16) tempPlus[(unsigned long)(3)];
            WorkPrec_ HONEI_ALIGNED(16) tempTopPlus[(unsigned long)(3)];

            WorkPrec_ HONEI_ALIGNED(16) phiPlusOld[(unsigned long)(3)];
            WorkPrec_ HONEI_ALIGNED(16) phiPlusNew[(unsigned long)(3)];
            WorkPrec_ HONEI_ALIGNED(16) phiMinusNew[(unsigned long)(3)];

            WorkPrec_ HONEI_ALIGNED(16) tempMinus[(unsigned long)(3)];
            WorkPrec_ HONEI_ALIGNED(16) tempTopMinus[(unsigned long)(3)];

            for(int i(0); i < 3; ++i)
            {
                tempPlus[i] = WorkPrec_(0);
                tempTopPlus[i] = WorkPrec_(0);
                phiPlusOld[i] = WorkPrec_(0);
                phiPlusNew[i] = WorkPrec_(0);
                phiMinusNew[i] = WorkPrec_(0);
                tempMinus[i] = WorkPrec_(0);
                tempTopMinus[i] = WorkPrec_(0);
            }

            WorkPrec_ HONEI_ALIGNED(16) phiMinusOld;
            WorkPrec_ HONEI_ALIGNED(16) temp;
            WorkPrec_ HONEI_ALIGNED(16) tempTop;
            WorkPrec_ HONEI_ALIGNED(16) prefac = delta_t/(4*delta_y);

            for(unsigned long s=0; s< d_width; ++s)
            {

                ///Needed Iterators.
                //typename DenseVector<WorkPrec_>::ElementIterator d(m2diag.begin_elements());
                unsigned long d(0);
                //typename DenseVector<WorkPrec_>::ConstElementIterator i(u.begin_elements());
                unsigned long i(0);
                //typename DenseVector<WorkPrec_>::ElementIterator b1(m2bandPlus1.begin_elements());
                unsigned long b1(0);
                //typename DenseVector<WorkPrec_>::ElementIterator b2(m2bandPlus2.begin_elements());
                unsigned long b2(0);
                //typename DenseVector<WorkPrec_>::ElementIterator bminus1(m2bandMinus1.begin_elements());
                unsigned long bminus1(0);

                ///Iterate to the next column
                for(unsigned long f=0; f < 3*(s+2); ++f)
                {
                    ++i;++d;++b1;++b2;++bminus1;
                }



                for(unsigned long k=0; k<3; ++k)
                {
                    tempPlus[k]= (w_e)[i] + (d_e)[k]*(u_e)[i];
                    ++i;++d;++b1;++b2;++bminus1;
                }

                //Iterate to next column-element
                for(unsigned long f=0; f<3*(d_width+3);++f)
                {
                    ++i;++d;++b1;++b2;++bminus1;
                }

                for(unsigned long k=0; k<3; ++k)
                {
                    temp= (w_e)[i] + (d_e)[k]*(u_e)[i];
                    tempTopPlus[k] = temp - tempPlus[k];
                    tempPlus[k] = temp;
                    tempMinus[k] =  (w_e)[i] - (d_e)[k]*(u_e)[i];
                    ++i;++d;++b1;++b2;++bminus1;
                }

                //Iterate i to next column-element
                for(unsigned long f=0; f<3*(d_width+3);++f)
                {
                    ++i;
                }

                for(unsigned long k=0; k<3; ++k)
                {
                    temp= (w_e)[i] - (d_e)[k]*(u_e)[i];
                    tempTopMinus[k] = temp - tempMinus[k];
                    tempMinus[k] = temp;
                    temp= (w_e)[i] + (d_e)[k]*(u_e)[i];
                    tempTop = temp - tempPlus[k];
                    if(tempTop != 0)
                    {
                        phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                    }
                    else
                    {
                        phiPlusOld[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                    }
                    tempPlus[k]=temp;
                    tempTopPlus[k]=tempTop;
                    ++i;
                }

                //Iterate i to next column-element
                for(unsigned long f=0; f<3*(d_width+3);++f)
                {
                    ++i;
                }

                for(unsigned long k=0; k<3; ++k)
                {
                    temp = (w_e)[i] - (d_e)[k]*((u_e)[i]);
                    tempTop = temp - tempMinus[k];
                    if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon())
                    {
                        phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);
                    }
                    else
                    {
                        phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                    }
                    tempMinus[k]=temp;
                    tempTopMinus[k] = tempTop;
                    temp  = (w_e)[i] + (d_e)[k]* (u_e)[i];
                    tempTop = temp - tempPlus[k];
                    if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon())
                    {
                        phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                    }
                    else
                    {
                        phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                    }
                    tempPlus[k]= temp;
                    tempTopPlus[k] = tempTop;
                    ++i;
                }
                for(unsigned long x = 0; x < d_height; ++x)
                {

                    //Iterate to next column-elements
                    for(unsigned long f=0; f<3*(d_width+3);++f)
                    {
                        ++i;++d;++b1;++b2;++bminus1;
                    }

                    for(unsigned long k =0; k<3; ++k)
                    {
                        temp = prefac *(2 - phiPlusOld[k]);

                        m2bandMinus1_e[bminus1] =temp;
                        m4bandMinus1_e[bminus1] =temp * (d_e)[k];
                        m2diag[d] = prefac * (phiPlusNew[k] + phiPlusOld[k] + phiMinusNew[k]);
                        m4diag[d] = (d_e)[k] * prefac *(phiPlusNew[k] + phiPlusOld[k] - phiMinusNew[k] - 4);

                        phiPlusOld[k]= phiPlusNew[k];
                        phiMinusOld = phiMinusNew[k];
                        temp = (w_e)[i] - (d_e)[k]*(u_e)[i];
                        tempTop = temp - tempMinus[k];
                        if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon())
                        {
                            phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/tempTop);
                        }
                        else
                        {
                            phiMinusNew[k] = min_mod_limiter(tempTopMinus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                        }
                        tempMinus[k]=temp;
                        tempTopMinus[k] = tempTop;
                        temp  = (w_e)[i] + (d_e)[k]* (u_e)[i];
                        tempTop = temp - tempPlus[k];
                        if( fabs(tempTop) >= std::numeric_limits<WorkPrec_>::epsilon())
                        {
                            phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/tempTop);
                        }
                        else
                        {
                            phiPlusNew[k] = min_mod_limiter(tempTopPlus[k]/std::numeric_limits<WorkPrec_>::epsilon());

                        }
                        tempPlus[k]= temp;
                        tempTopPlus[k] = tempTop;
                        m2bandPlus1_e[b1] = prefac * (-2 - phiPlusOld[k] - phiMinusOld - phiMinusNew[k]);
                        m4bandPlus1_e[b1] = (d_e)[k]*prefac * (2 - phiPlusOld[k] + phiMinusOld + phiMinusNew[k]);
                        m2bandPlus2_e[b2] = prefac* phiMinusNew[k];
                        m4bandPlus2_e[b2] = (d_e)[k] * prefac *(-phiMinusNew[k]);
                        ++i;++d;++b1;++b2;++bminus1;
                    }
                }
            }
            m2.insert_band(0, m2diag);
            m2.insert_band(3*(d_width+4), m2bandPlus1);
            m2.insert_band(6*(d_width+4), m2bandPlus2);
            m2.insert_band((-3)*(d_width+4), m2bandMinus1);
            m4.insert_band(0, m4diag);
            m4.insert_band(3*(d_width+4), m4bandPlus1);
            m4.insert_band(6*(d_width+4), m4bandPlus2);
            m4.insert_band((-3)*(d_width+4), m4bandMinus1);
#ifdef SOLVER_VERBOSE
            std::cout << "Finished Matrix Assembly 2.\n";
#endif
            }
    };

    template <typename Tag_>
        struct AssemblyProcessing<Tag_, assembly_types::QUICK::M6>
        {

            public:
                template <typename WorkPrec_>
                    static inline BandedMatrix<WorkPrec_> value(BandedMatrix<WorkPrec_> & m1, BandedMatrix<WorkPrec_> & result, DenseVector<WorkPrec_> & c, unsigned long d_width, unsigned long d_height)
                    {
                        CONTEXT("When processing RelaxSolver assembly, quick stage, matrix m_6:");

                        ///Bands of the matrix which will be assembled.
                        DenseVector<WorkPrec_> m6diag = (m1.band((unsigned long)(0))).copy();
                        DenseVector<WorkPrec_> m6bandplus3 = (m1.band((unsigned long)(3))).copy();
                        DenseVector<WorkPrec_> m6bandplus6 = (m1.band((unsigned long)(6))).copy();
                        DenseVector<WorkPrec_> m6bandminus3 = (m1.band((unsigned long)(-3))).copy();
                        ///Needed Iterators.
                        typename DenseVector<WorkPrec_>::ElementIterator d(m6diag.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator b1(m6bandplus3.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator b2(m6bandplus6.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator bminus1(m6bandminus3.begin_elements());

                        DenseVector<WorkPrec_> c_squared((c.copy()));
                        ElementProduct<Tag_>::value(c_squared, (c));

                        for( ; d.index() < 6*(d_width+4); ++d);
                        for( ; b1.index() < 6*(d_width+4); ++b1);
                        for( ; b2.index() < 6*(d_width+4); ++b2);
                        for( ; bminus1.index() < 6*(d_width+4); ++bminus1);

                        while(d.index() < 3*(d_width+4)*(d_height+2))
                        {
                            ++d; ++d; ++d; ++d; ++d; ++d;
                            ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
                            ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
                            ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;

                            for(unsigned long i = 0; i < d_width; ++i)
                            {
                                *d *= c_squared[0];
                                *b1 *= c_squared[0];
                                *b2 *= c_squared[0];
                                *bminus1 *= c_squared[0];
                                ++d; ++b1; ++b2; ++bminus1;
                                *d *= c_squared[1];
                                *b1 *= c_squared[1];
                                *b2 *= c_squared[1];
                                *bminus1 *= c_squared[1];
                                ++d; ++b1; ++b2; ++bminus1;
                                *d *= c_squared[2];
                                *b1 *= c_squared[2];
                                *b2 *= c_squared[2];
                                *bminus1 *= c_squared[2];
                                ++d; ++b1; ++b2; ++bminus1;
                            }

                            ++d; ++d; ++d; ++d; ++d; ++d;
                            ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
                            ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
                            ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;
                        }
                        result.insert_band(0, m6diag);
                        result.insert_band(3, m6bandplus3);
                        result.insert_band(6, m6bandplus6);
                        result.insert_band(-3, m6bandminus3);
#ifdef SOLVER_VERBOSE
                        std::cout << "Finished Quick Assembly m2.\n";
#endif
                        return result;

                    }
        };
    template <typename Tag_>
        struct AssemblyProcessing<Tag_, assembly_types::QUICK::M8>
        {

            public:
                template <typename WorkPrec_>
                    static inline BandedMatrix<WorkPrec_> value(BandedMatrix<WorkPrec_> & m2, BandedMatrix<WorkPrec_> & result, DenseVector<WorkPrec_> & dv, unsigned long d_width, unsigned long d_height)
                    {
                        CONTEXT("When processing RelaxSolver assembly, quick stage, matrix m_8:");

                        DenseVector<WorkPrec_> m8diag = (m2.band((unsigned long)(0))).copy();
                        DenseVector<WorkPrec_> m8bandplus3 = (m2.band((unsigned long)(3*(d_width +4)))).copy();
                        DenseVector<WorkPrec_> m8bandplus6 = (m2.band((unsigned long)(6*(d_width +4)))).copy();
                        DenseVector<WorkPrec_> m8bandminus3 = (m2.band((unsigned long)((-3)*(d_width +4)))).copy();
                        ///Needed Iterators.
                        typename DenseVector<WorkPrec_>::ElementIterator d(m8diag.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator b1(m8bandplus3.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator b2(m8bandplus6.begin_elements());
                        typename DenseVector<WorkPrec_>::ElementIterator bminus1(m8bandminus3.begin_elements());
                        DenseVector<WorkPrec_> d_squared(((dv.copy())));
                        ElementProduct<Tag_>::value(d_squared, (dv));

                        for( ; d.index() < 6*(d_width + 4); ++d);
                        for( ; b1.index() < 6*(d_width + 4); ++b1);
                        for( ; b2.index() < 6*(d_width + 4); ++b2);
                        for( ; bminus1.index() < 6*(d_width + 4); ++bminus1);

                        while(d.index() < 3*(d_width+4)*(d_height+2))
                        {
                            ++d; ++d; ++d; ++d; ++d; ++d;
                            ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
                            ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
                            ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;

                            for(unsigned long i = 0; i < d_width; ++i)
                            {
                                *d *= d_squared[0];
                                *b1 *= d_squared[0];
                                *b2 *= d_squared[0];
                                *bminus1 *= d_squared[0];
                                ++d; ++b1; ++b2; ++bminus1;
                                *d *= d_squared[1];
                                *b1 *= d_squared[1];
                                *b2 *= d_squared[1];
                                *bminus1 *= d_squared[1];
                                ++d; ++b1; ++b2; ++bminus1;
                                *d *= d_squared[2];
                                *b1 *= d_squared[2];
                                *b2 *= d_squared[2];
                                *bminus1 *= d_squared[2];
                                ++d; ++b1; ++b2; ++bminus1;
                            }

                            ++d; ++d; ++d; ++d; ++d; ++d;
                            ++b1; ++b1; ++b1; ++b1; ++b1; ++b1;
                            ++b2; ++b2; ++b2; ++b2; ++b2; ++b2;
                            ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1; ++bminus1;
                        }
                        result.insert_band(0, m8diag);
                        result.insert_band(3*(d_width +4), m8bandplus3);
                        result.insert_band(6*(d_width +4), m8bandplus6);
                        result.insert_band((-3)*(d_width +4),m8bandminus3);
#ifdef SOLVER_VERBOSE
                        std::cout << "Finished Quick Assembly m4.\n";
#endif
                        return result;
                    }
        };

    //------------------SSE--------------------
    template<>
        struct AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M6>
        {
            public:
                static BandedMatrix<float> & value(BandedMatrix<float> & m1, BandedMatrix<float> & result, DenseVector<float> & c, unsigned long d_width, unsigned long d_height);
                static BandedMatrix<double> & value(BandedMatrix<double> & m1, BandedMatrix<double> & result, DenseVector<double> & c, unsigned long d_width, unsigned long d_height);
        };

    template<>
        struct AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M8>
        {
            public:
                static BandedMatrix<float> & value(BandedMatrix<float> & m2, BandedMatrix<float> & result, DenseVector<float> & c, unsigned long d_width, unsigned long d_height);
                static BandedMatrix<double> & value(BandedMatrix<double> & m2, BandedMatrix<double> & result, DenseVector<double> & c, unsigned long d_width, unsigned long d_height);
        };
//--------------------------MCSSE------------------------------------------
    template<>
        struct AssemblyProcessing<tags::CPU::MultiCore::SSE, assembly_types::QUICK::M6>
        {
            public:
                static BandedMatrix<float> & value(BandedMatrix<float> & m1, BandedMatrix<float> & result, DenseVector<float> & c, unsigned long d_width, unsigned long d_height)
                {
                    return AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M6>::value(m1, result, c, d_width, d_height);
                }
                static BandedMatrix<double> & value(BandedMatrix<double> & m1, BandedMatrix<double> & result, DenseVector<double> & c, unsigned long d_width, unsigned long d_height)
                {
                    return AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M6>::value(m1, result, c, d_width, d_height);

                }
        };

    template<>
        struct AssemblyProcessing<tags::CPU::MultiCore::SSE, assembly_types::QUICK::M8>
        {
            public:
                static BandedMatrix<float> & value(BandedMatrix<float> & m2, BandedMatrix<float> & result, DenseVector<float> & c, unsigned long d_width, unsigned long d_height)
                {
                    return AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M8>::value(m2, result, c, d_width, d_height);

                }
                static BandedMatrix<double> & value(BandedMatrix<double> & m2, BandedMatrix<double> & result, DenseVector<double> & c, unsigned long d_width, unsigned long d_height)
                {
                    return AssemblyProcessing<tags::CPU::SSE, assembly_types::QUICK::M8>::value(m2, result, c, d_width, d_height);
                }
        };

}
#endif
