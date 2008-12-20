/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de.de>
 *
 * This file is part of the LibVisual C++ library. LibVisual is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibVisual is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBSWE_GUARD_SOLVER_SERVER__HH
#define LIBSWE_GUARD_SOLVER_SERVER_HH 1
#include <honei/la/dense_matrix.hh>
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <honei/la/dense_matrix.hh>
#include <honei/swe/relax_solver.hh>
#include <honei/swe/volume.hh>
#include <cstring>

#define BUFFER_SIZE 1024

namespace honei
{
    namespace globals {

        //globally defined solver:
        ulint dwidth = 41;
        ulint dheight = 41;
        //DenseMatrix<float> height(dheight, dwidth, float(5));
        DenseMatrix<float> height(dheight, dwidth, float(5));

        DenseMatrix<float> bottom(dheight, dwidth, float(0));
        DenseMatrix<float> u1(dheight, dwidth, float(0));
        DenseMatrix<float> u2(dheight, dwidth, float(0));
        unsigned long entries = 3*((dwidth*dheight)+4*(dwidth+dheight+4));
        DenseVector<float> u(entries, float(1));
        DenseVector<float> v(entries, float(1));
        DenseVector<float> w(entries, float(1));
        DenseVector<float> bx (entries/3, float(1));
        DenseVector<float> by (entries/3, float(1));
        DenseVector<float> c (3,float(5));
        DenseVector<float> d (3,float(5));
        float deltax = 5;
        float deltay = 5;
        float deltat = 5./24.;
        double eps = 10e-6;
        float manning = float(0);
#ifdef HONEI_SSE
        RelaxSolver<tags::CPU::SSE, float, float, float, float, float, SIMPLE, REFLECT, precision_modes::FIXED> solver( &height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);
#elif defined HONEI_CELL
        RelaxSolver<tags::Cell, float, float, float, float, float, SIMPLE, REFLECT, precision_modes::FIXED> solver( &height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);
#else
        RelaxSolver<tags::CPU, float, float, float, float, float, SIMPLE, REFLECT, precision_modes::FIXED> solver( &height, &bottom, &u1, &u2, &u, &v, &w,
                dwidth, dheight, deltax, deltay, deltat, eps, &bx, &by, &c, &d, manning);
#endif

    }

    template <typename Tag_, typename DataType_> class SolverServer
    {
        private:
            int _scenario;

            void _read_scenario(int c)
            {
                //char buffer[BUFFER_SIZE];
                char scenario[BUFFER_SIZE];
                int bytes;

                bytes = recv(c, scenario, sizeof(scenario) - 1, 0);
                scenario[bytes] = '\0';
                _scenario = atoi(scenario);
                //std::cout<<"scenario choosen: "<<_scenario<<std::endl;
            }

            int _write_timesteps(int c)
            {
                char handshake[BUFFER_SIZE] ;
                int bytes;

                std::string result;
                //initial scenario setup
                Cylinder<float> c1(globals::height, float(15.), globals::dwidth/2, globals::dheight/2);
                c1.value();

                globals::c[0] = 12;
                globals::c[1] = 7;
                globals::c[2] = 12;
                globals::d[0] = 12;
                globals::d[1] = 7;
                globals::d[2] = 12;
                globals::solver.do_preprocessing();

                do
                {
                    //std::cout<<"Timestep:"<<std::endl;
                    // insert data calculation by solver here
                    globals::solver.solve();
                    for (unsigned long row(0) ; row < globals::height.rows() ; ++row)
                    {
                        result = "";
                        for (typename DenseVector<DataType_>::ElementIterator i((globals::height)[row].begin_elements()),
                                i_end((globals::height)[row].end_elements()) ; i != i_end ; ++i)
                        {
                            result += stringify(*i) + "a#";
                        }
                        //std::cout<< "sending matrix row: "<<(globals::height)[row]<<std::endl;
                        bytes = send(c, result.c_str(), strlen(result.c_str()) - 2, 0);

                        bytes = recv(c, handshake, sizeof(handshake) - 1, 0);
                        handshake[bytes] = '\0';
                        if (handshake[0] == 'r' || handshake[0] == 'q' || handshake[0] == 'x' || handshake[0] == 'e')
                            break;
                    }
                }
                while (handshake[0] != 'r' && handshake[0] != 'q' && handshake[0] != 'x');
                if (handshake[0] == 'r')
                    return 1;
                else if (handshake[0] == 'q')
                    return 0;
                else if (handshake[0] == 'x')
                    return 2;
                else
                    throw InternalError("Invalid handhake status!");
            }

            int _handling(int c)
            {
                _read_scenario(c);
                int retval(0);
                do
                {
                    std::cout<<"Starting timesteps..."<<std::endl;
                    retval = _write_timesteps(c);
                }
                while (retval == 1);
                return retval;
            }

        public:
            SolverServer()
            {
            }

            ~SolverServer()
            {
            }

            void run()
            {
                int s, c, cli_size;
                struct sockaddr_in srv, cli;

                s = socket(AF_INET, SOCK_STREAM, 0);
                if (s == -1)
                {
                    perror("socket() failed");
                }

                srv.sin_addr.s_addr = INADDR_ANY;
                srv.sin_port = htons(4711);
                srv.sin_family = AF_INET;

                if (bind(s, (struct sockaddr*)&srv, sizeof(srv)) == -1)
                {
                    perror("bind() failed");
                }

                if (listen(s, 3) == -1)
                {
                    perror("listen() failed");
                }

                for(;;)
                {
                    std::cout<<"Waiting for clients to connect."<<std::endl;
                    cli_size = sizeof(cli);
                    c = accept(s, (struct sockaddr*)&cli, (socklen_t*)&cli_size);
                    if (c == -1)
                    {
                        perror("accept() failed");
                    }

                    if (_handling(c) == 2)
                    {
                        //shutdown server
                        close(c);
                        return;
                    }
                    else
                        close(c);
                }
            }

    };

}
#endif
