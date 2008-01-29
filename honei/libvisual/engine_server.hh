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
#include <honei/libla/dense_matrix.hh>
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
#include <honei/libla/dense_matrix.hh>

#define BUFFER_SIZE 1024

namespace honei
{
    template <typename Tag_, typename DataType_> class EngineServer
    {
        private:
            int _scenario;
            DenseMatrix<DataType_> * _height_field;

            void _read_scenario(int c)
            {
                char buffer[BUFFER_SIZE], scenario[BUFFER_SIZE];
                int bytes;

                bytes = recv(c, scenario, sizeof(scenario) - 1, 0);
                scenario[bytes] = '\0';
                _scenario = atoi(scenario);

                delete _height_field;
                _height_field = new DenseMatrix<DataType_>(4, 4, DataType_(0));
                for (typename MutableMatrix<DataType_>::ElementIterator i(_height_field->begin_elements()),
                        i_end(_height_field->end_elements()) ; i != i_end ; ++i)
                {
                    *i = DataType_(i.index() + 1) / DataType_(1.123);
                }
                std::cout<<"scenario choosen: "<<_scenario<<std::endl;
            }

            int _write_timesteps(int c)
            {
                char handshake[BUFFER_SIZE] ;
                int bytes;

                std::string result;

                do
                {
                    std::cout<<"Timestep:"<<std::endl;
                    // insert data calculation by solver here
                    for (unsigned long row(0) ; row < _height_field->rows() ; ++row)
                    {
                        result = "";
                        for (typename DenseVector<DataType_>::ElementIterator i((*_height_field)[row].begin_elements()),
                                i_end((*_height_field)[row].end_elements()) ; i != i_end ; ++i)
                        {
                            *i = *i + 1.1 * row;
                            result += stringify(*i) + "a#";
                        }
                        std::cout<< "sending matrix row: "<<(*_height_field)[row]<<std::endl;
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
            EngineServer() :
                _height_field(0)
            {
            }

            ~EngineServer()
            {
                delete _height_field;
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
