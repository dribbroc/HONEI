/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de.de>
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

#pragma once
#ifndef LIBSWE_GUARD_SOLVER_CLIENT_HH
#define LIBSWE_GUARD_SOLVER_CLIENT_HH 1
#include <honei/la/dense_matrix.hh>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <honei/la/dense_matrix.hh>
#include <cstring>
#include <string>

#define BUFFER_SIZE 1024
#define DATA_SIZE 5234567

namespace honei
{
    template <typename Tag_, typename DataType_> class SolverClient
    {
        private:
            int _socket;
            char _data[DATA_SIZE];

            void _write_scenario(int c, int scenario)
            {
                char buffer[BUFFER_SIZE];
                int bytes;


                //std::cout<<"Selecting scenario"<<std::endl;
                sprintf(buffer, "%i", scenario);
                bytes = send(c, buffer, strlen(buffer), 0);
            }

            void _read_timestep(int c, DenseMatrix<double> & height_field)
            {
                int bytes;
                char delims[] = "#";
                typename DenseMatrix<double>::ElementIterator i(height_field.begin_elements()), i_end(height_field.end_elements());

                for (unsigned long row(0) ; row < height_field.rows() ; ++row)
                {
                    bytes = recv(c, _data, sizeof(_data) - 1, 0);
                    _data[bytes] = '\0';

                    char *result = NULL;
                    result = strtok(_data, delims);
                    while (result != NULL && i != i_end)
                    {
                        *i = atof(result);
                        result = strtok(NULL, delims);
                        ++i;
                    }

                    bytes = send(c, "c", 1, 0);
                }
                //std::cout<<"got matrix: "<<height_field<<std::endl;
            }

            void _restart(int c)
            {
                int bytes;
                bytes = recv(c, _data, sizeof(_data) - 1, 0);
                //std::cout<<"Restarting scenario."<<std::endl;
                send(c, "r", 1 ,0);
            }

            void _quit(int c)
            {
                send(c, "q", 1 ,0);
            }

            void _shutdown(int c)
            {
                send(c, "x", 1, 0);
            }

        public:
            SolverClient() :
                _socket(-1)
            {
            }

            ~SolverClient()
            {
                if (_socket != -1) close(_socket);
            }

            void init(const char *  hostname, int port, int scenario)
            {
                if (_socket != -1) close(_socket);
                struct sockaddr_in srv;

                _socket = socket(AF_INET, SOCK_STREAM, 0);
                if (_socket == -1)
                {
                    perror("socket failed()");
                }

                srv.sin_addr = *(struct in_addr*) gethostbyname(hostname)->h_addr;

                srv.sin_port = htons(port);
                srv.sin_family = AF_INET;

                if (connect(_socket, (struct sockaddr*)&srv, sizeof(srv)) == -1)
                {
                    perror("connect failed()");
                }

                _write_scenario(_socket, 12345);

            }

            DenseMatrix<double> & do_step (DenseMatrix<double> & height_field)
            {
                _read_timestep(_socket, height_field);
                return height_field;
            }

            void restart_scenario()
            {
                _restart(_socket);
            }

            void quit_server()
            {
                _quit(_socket);
                if (_socket != -1) close(_socket);
            }

            void shutdown_server()
            {
                _shutdown(_socket);
            }
    };
}
#endif
