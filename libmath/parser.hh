/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBUTIL_GUARD_PARSER_HH
#define LIBUTIL_GUARD_PARSER_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

namespace honei
{
    /**
     * \brief Parser for DenseVectors and Matrices.
     *
     * The fileformat is straightforward: DenseVectors are stored serially in one
     * record, DenseMatrices are stores rowwisely, each record separated by a single
     * '/' character.
     *
     * \ingroup grpdebug
     */


    template<typename DT_>
    struct Parser
    {

    };

    template<>
    struct Parser<float>
    {
        public:
            /**
            * \brief Returns floating point number from string specified.
            *
            * \param c The C String to convert.
            *
            */

            /// \{
            static float parse(const char* c)
            {
                CONTEXT("When parsing string.");
                return float(std::atof(c));
            }

            /**
            * \brief Returns a DenseVector read and parsed from file.
            *
            * \param filename The name of the file to read in.
            * \param size The length of the vector.
            *
            */

            /// \{
            static DenseVector<float> parse(char* filename, unsigned long size)
            {
                CONTEXT("When parsing file to DenseVector");
                DenseVector<float> result(size, float(0.0));
                ifstream ifs;
                ifs.open(filename);
                char * line = new char[size*100];

                for(long a = 0; a < size*100; a++)
                {
                    line[a] = ' ';
                }

                ifs.getline(line, size*100);

                DenseVector<float>::ElementIterator i(result.begin_elements()), i_end(result.end_elements());

                char * target = new char [size*100];
                for(long a = 0; a < size*100; a++)
                {
                    target[a] = ' ';
                }

                for(unsigned long j = 0; j < size *100; j++)
                {
                    if(line[j] != ' ' && ifs.good())
                    {
                        target[j] = line[j];
                    }
                    else
                        if(i != i_end)
                        {
                            *i = parse(target);
                            target = new char [size*100];
                            for(long a = 0; a < size*100; a++)
                            {
                                target[a] = ' ';
                            }


                            ++i;
                        }
                }
                delete target;
                delete line;

                ifs.close();
                return result;

            }

            /**
            * \brief Returns a DenseMatrix read and parsed from file.
            *
            * \param filename The name of the file to read in.
            * \param width The width of the matrix.
            * \param filename The height of the matrix.
            *
            */

            /// \{
            static DenseMatrix<float> parse(char* filename, unsigned long width, unsigned long height)
            {
                CONTEXT("When parsing file to DenseMatrix");
                DenseMatrix<float> result(height, width, float(0.0));
                char lines[height][width*100];
                unsigned long linescount = 0;
                ifstream ifs;
                ifs.open(filename);
                string s;
                while(getline(ifs, s))
                {
                    for(unsigned long i = 0; i < width*100; i++)
                    {
                        lines[linescount][i] = s[i];
                    }
                    linescount++;
                }
                ifs.close();
                for(unsigned long i = 0; i < height; ++i)
                {
                    unsigned long lineprogress = 0;
                    unsigned int tempprogress = 0;

                    string s1;
                    for(unsigned long j = 0; j < width*100; ++j)
                    {
                        if(lines[i][j] == ' ' || lines[i][j] == '/')
                        {
                            if(lineprogress < width) 
                            {
                                result[i][lineprogress] = Parser<float>::parse(s1.c_str());
                                ++lineprogress;
                                tempprogress = 0;
                                for(int a = 0; a < width*100; a++)
                                {
                                    s1[a] = ' ';
                                }
                            }

                        }
                        else
                        {
                            s1[tempprogress] = lines[i][j];
                            ++tempprogress;
                        }
                    }
                }
                return result;
            }

    };

    template<>
    struct Parser<double>
    {
        public:
            /**
            * \brief Returns floating point number from string specified.
            *
            * \param c The C String to convert.
            *
            */

            /// \{


            static double parse(const char* c)
            {
                CONTEXT("When parsing string.");
                return double(std::atof(c));
            }
            /**
            * \brief Returns a DenseVector read and parsed from file.
            *
            * \param filename The name of the file to read in.
            * \param size The length of the vector.
            *
            */

            /// \{
            static DenseVector<double> parse(char* filename, unsigned long size)
            {
                CONTEXT("When parsing file to DenseVector");
                DenseVector<double> result(size, double(0.0));
                ifstream ifs;
                ifs.open(filename);
                char * line = new char[size*100];

                for(long a = 0; a < size*100; a++)
                {
                    line[a] = ' ';
                }

                ifs.getline(line, size*100);

                DenseVector<double>::ElementIterator i(result.begin_elements()), i_end(result.end_elements());

                char * target = new char [size*100];
                for(long a = 0; a < size*100; a++)
                {
                    target[a] = ' ';
                }

                for(unsigned long j = 0; j < size *100; j++)
                {
                    if(line[j] != ' ' && ifs.good())
                    {
                        target[j] = line[j];
                    }
                    else
                        if(i != i_end)
                        {
                            *i = parse(target);
                            target = new char [size*100];
                            for(long a = 0; a < size*100; a++)
                            {
                                target[a] = ' ';
                            }


                            ++i;
                        }
                }
                delete target;
                delete line;

                ifs.close();
                return result;

            }
            /**
            * \brief Returns a DenseMatrix read and parsed from file.
            *
            * \param filename The name of the file to read in.
            * \param width The width of the matrix.
            * \param filename The height of the matrix.
            *
            */

            /// \{
            static DenseMatrix<double> parse(char* filename, unsigned long width, unsigned long height)
            {
                CONTEXT("When parsing file to DenseMatrix");
                DenseMatrix<double> result(height, width, double(0.0));
                char lines[height][width*50];
                unsigned long linescount = 0;
                ifstream ifs;
                ifs.open(filename);
                string s;
                while(getline(ifs, s))
                {
                    for(unsigned long i = 0; i < width*50; i++)
                    {
                        lines[linescount][i] = s[i];
                    }
                    linescount++;
                }
                ifs.close();

                for(unsigned long i = 0; i < height; ++i)
                {
                    unsigned long lineprogress = 0;
                    unsigned int tempprogress = 0;

                    string s1;
                    for(unsigned long j = 0; j < width*50; ++j)
                    {
                        if(lines[i][j] == ' ' || lines[i][j] == '/')
                        {
                            if(lineprogress < width)
                            {
                                result(i,lineprogress) = Parser<double>::parse(s1.c_str());
                                ++lineprogress;
                                tempprogress = 0;
                                for(int a = 0; a < width*50; a++)
                                {
                                    s1[a] = ' ';
                                }
                            }

                        }
                        else
                        {
                            s1[tempprogress] = lines[i][j];
                            ++tempprogress;
                        }
                    }
                }
                return result;
            }

    };

}
#endif
