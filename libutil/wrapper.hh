/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Joachim Messer <joachim.messer@uni-dortmund.de>
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_WRAPPER_HH
#define LIBLA_GUARD_WRAPPER_HH 1

#include <libutil/lock.hh>

namespace honei
{
    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_>
        class ResultThreeArgWrapper
        {
            private:
                DT1_ & _result;

                DT2_ _a;

                DT3_ _b;

                DT4_ _c;


            public:

                typedef void result_type;

                ResultThreeArgWrapper(DT1_ & result, DT2_ & a, DT3_ & b, DT4_ & c):
                    _result(result),
                    _a(a),
                    _b(b),
                    _c(c)
                {
                }

                void operator() ()
                {
                    _result = Tag_::value(_a, _b, _c);
                }
                
                void operator()(Mutex * mutex)
                {
                    DT1_ temp(Tag_::value(_a, _b, _c));
                    Lock l(*mutex);
                    _result += temp;
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_>
        class ThreeArgWrapper
        {
            private:
                DT1_ _a;

                DT2_ _b;

                DT3_ _c;


            public:

                typedef void result_type;

                ThreeArgWrapper(DT1_ & a, DT2_ & b, DT3_ & c):
                    _a(a),
                    _b(b),
                    _c(c)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a, _b, _c);
                }

        };
    
    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_>
        class ResultTwoArgWrapper
        {
            private:
                DT1_ & _result;

                DT2_ _a;

                DT3_ _b;


            public:

                typedef void result_type;

                ResultTwoArgWrapper(DT1_ & result, DT2_ & a, DT3_ & b):
                    _result(result),
                    _a(a),
                    _b(b)
                {
                }

                void operator() ()
                {
                    _result = Tag_::value(_a, _b);
                }
                
                void operator()(Mutex * mutex)
                {
                    DT1_ temp(Tag_::value(_a, _b));
                    Lock l(*mutex);
                    _result += temp;
                }

        };


    template < typename Tag_, typename DT2_, typename DT3_>
        class TwoArgWrapper
        {
            private:
                DT2_ _a;

                DT3_ _b;


            public:

                typedef void result_type;

                TwoArgWrapper(DT2_ & a, DT3_ & b):
                    _a(a),
                    _b(b)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a, _b);
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_>
        class ResultOneArgWrapper
        {
            private:
                DT1_ & _result;

                DT2_ _a;

            public:

                typedef void result_type;

                ResultOneArgWrapper(DT1_ & result, DT2_ & a):
                    _result(result),
                    _a(a)
                {
                }

                void operator() ()
                {
                    _result = Tag_::value(_a);
                }
                
                void operator()(Mutex * mutex)
                {
                    DT1_ temp(Tag_::value(_a));
                    Lock l(*mutex);
                    _result += *temp;
                }

        };

    template < typename Tag_, typename DT1_>
        class OneArgWrapper
        {
            private:
                DT1_ & _a;

            public:

                typedef void result_type;

                OneArgWrapper(DT1_ & a):
                    _a(a)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a);
                }

        };


}
#endif


