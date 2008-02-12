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

#include <honei/libutil/lock.hh>
#include <honei/libutil/profiler.hh>

namespace honei
{
    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_, typename DT5_, typename DT6_>
        class SixArgWrapper
        {
            private:
                DT1_ _a;

                DT2_ _b;

                DT3_ _c;

                DT4_ _d;

                DT5_ _e;

                DT6_ _f;


            public:

                typedef void result_type;

                SixArgWrapper(DT1_ & a, DT2_ & b, DT3_ & c, DT4_ & d, DT5_ & e, DT6_ & f):
                    _a(a),
                    _b(b),
                    _c(c),
                    _d(d),
                    _e(e),
                    _f(f)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a, _b, _c, _d, _e, _f);
                }

                void operator()(Mutex * mutex)
                {
                    PROFILER_START("Wrapper<6>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<6>:mutex");
                    Tag_::value(_a, _b, _c, _d, _e, _f);
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_, typename DT5_, typename DT6_>
        class ResultFiveArgWrapper
        {
            private:
                DT1_ & _result;

                DT2_ _a;

                DT3_ _b;

                DT4_ _c;

                DT5_ _d;

                DT6_ _e;


            public:

                typedef void result_type;

                ResultFiveArgWrapper(DT1_ & result, DT2_ & a, DT3_ & b, DT4_ & c, DT5_ & d, DT6_ & e):
                    _result(result),
                    _a(a),
                    _b(b),
                    _c(c),
                    _d(d),
                    _e(e)
                {
                }

                void operator() ()
                {
                    _result = Tag_::value(_a, _b, _c, _d, _e);
                }

                void operator()(Mutex * mutex)
                {
                    DT1_ temp(Tag_::value(_a, _b, _c, _d, _e));
                    PROFILER_START("Wrapper<5,R>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<5,R>:mutex");
                    _result += temp;
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_, typename DT5_>
        class FiveArgWrapper
        {
            private:
                DT1_ _a;

                DT2_ _b;

                DT3_ _c;

                DT4_ _d;

                DT5_ _e;


            public:

                typedef void result_type;

                FiveArgWrapper(DT1_ & a, DT2_ & b, DT3_ & c, DT4_ & d, DT5_ & e):
                    _a(a),
                    _b(b),
                    _c(c),
                    _d(d),
                    _e(e)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a, _b, _c, _d, _e);
                }

                void operator()(Mutex * mutex)
                {
                    PROFILER_START("Wrapper<5>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<5>:mutex");
                    Tag_::value(_a, _b, _c, _d, _e);
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_, typename DT5_>
        class ResultFourArgWrapper
        {
            private:
                DT1_ & _result;

                DT2_ _a;

                DT3_ _b;

                DT4_ _c;

                DT5_ _d;


            public:

                typedef void result_type;

                ResultFourArgWrapper(DT1_ & result, DT2_ & a, DT3_ & b, DT4_ & c, DT5_ & d):
                    _result(result),
                    _a(a),
                    _b(b),
                    _c(c),
                    _d(d)
                {
                }

                void operator() ()
                {
                    _result = Tag_::value(_a, _b, _c, _d);
                }

                void operator()(Mutex * mutex)
                {
                    DT1_ temp(Tag_::value(_a, _b, _c, _d));
                    PROFILER_START("Wrapper<4,R>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<4,R>:mutex");
                    _result += temp;
                }

        };

    template < typename Tag_, typename DT1_, typename DT2_, typename DT3_, typename DT4_>
        class FourArgWrapper
        {
            private:
                DT1_ _a;

                DT2_ _b;

                DT3_ _c;

                DT4_ _d;

            public:

                typedef void result_type;

                FourArgWrapper(DT1_ & a, DT2_ & b, DT3_ & c, DT4_ & d):
                    _a(a),
                    _b(b),
                    _c(c),
                    _d(d)
                {
                }

                void operator() ()
                {
                    Tag_::value(_a, _b, _c, _d);
                }

                void operator() (Mutex * mutex)
                {
                    PROFILER_START("Wrapper<4>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<4>:mutex");
                    Tag_::value(_a, _b, _c, _d);
                }
        };

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
                    PROFILER_START("Wrapper<3,R>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<3,R>:mutex");
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

                void operator() (Mutex * mutex)
                {
                    PROFILER_START("Wrapper<3>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<3>:mutex");
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
                    PROFILER_START("Wrapper<2,R>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<2,R>:mutex");
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
                    PROFILER_START("Wrapper<1,R>:mutex");
                    Lock l(*mutex);
                    PROFILER_STOP("Wrapper<1,R>:mutex");
                    _result += temp;
                }

        };

    template < typename Tag_, typename DT1_>
        class OneArgWrapper
        {
            private:
                DT1_ _a;

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


