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

#ifndef LIBSWE_GUARD_SOLVER_CLIENT__HH
#define LIBSWE_GUARD_SOLVER_CLIENT_HH 1
#include <honei/libla/dense_matrix.hh>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <honei/libla/dense_matrix.hh>

#define BUFFER_SIZE 1024
#define DATA_SIZE 1234567

namespace honei
{
    template <typename Tag_, typename DataType_> class SolverClient
    {
        private:
            int _socket;

            void _write_scenario(int c, int scenario)
            {
                char buffer[BUFFER_SIZE];
                int bytes;


                std::cout<<"Selecting scenario"<<std::endl;
                //strcpy(buffer, "4711");
                //itoa(scenario, buffer, 10);
                sprintf(buffer, "%i", scenario);
                bytes = send(c, buffer, strlen(buffer), 0);
            }

            void _read_timesteps(int c)
            {
                char buffer[DATA_SIZE];
                int bytes;

                unsigned count(0);
                do
                {
                    bytes = recv(c, buffer, sizeof(buffer) - 1, 0);
                    buffer[bytes] = '\0';
                    std::cout<<"Received tve: "<<buffer<<std::endl;

                    strcpy(buffer, "c");
                    bytes = send(c, buffer, 1, 0);

                    count++;
                }
                while (count < 100);
                strcpy(buffer, "q");
                bytes = send(c, buffer, 1, 0);
            }

            int _handling(int c)
            {
                //_write_scenario(c, 12345);
                _read_timesteps(c);
                return 0;
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

            void init(char *  hostname, int port, int scenario)
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
            void run()
            {
                _handling(_socket);
            }

            DenseMatrix<DataType_> & do_step (DenseMatrix<DataType_> & height_field)
            {
                return height_field;
            }
    };
}
#endif
