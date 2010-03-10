/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MULTICORE_GUARD_NUMAINFO_HH
#define MULTICORE_GUARD_NUMAINFO_HH 1

#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sched.h>
#include <dirent.h>

#define howmany(x,y) (((x)+((y)-1))/(y))
#define bitsperlong (8 * sizeof(unsigned long))
#define longsperbits(n) howmany(n, bitsperlong)
#define round_up(x,y) (((x) + (y) - 1) & ~((y)-1))

#define CPU_BYTES(x) (round_up(x, sizeof(unsigned long)))
#define CPU_LONGS(x) (CPU_BYTES(x) / sizeof(unsigned long))

#if defined(__x86_64__)
#define __NR_sched_getaffinity    204

#elif defined (__i386__)
#define __NR_sched_getaffinity    242

#elif defined(__ia64__)
#define __NR_sched_getaffinity    1232

#elif defined(__ppc__) || defined (__PPU__)
#define __NR_sched_getaffinity    233
#endif

namespace honei
{
    namespace mc
    {
        namespace intern
        {
            struct NumaInfo
            {
                struct bitmask
                {
                    unsigned long size;
                    unsigned long * maskp;
                };

                unsigned mask_size;
                unsigned num_nodes;
                struct bitmask ** node_cpu_mask;

                bool init(unsigned nodes)
                {
                    num_nodes = nodes;
                    cpu_set_t mask;
                    mask_size = 8 * syscall(__NR_sched_getaffinity, 0, sizeof(cpu_set_t), &mask);
                    node_cpu_mask = (bitmask **) calloc(num_nodes, sizeof(struct bitmask *));

                    return create_node_cpu_mask();
                }

                void cleanup()
                {
                    for (unsigned node(0) ; node < num_nodes ; ++node)
                    {
                        free(node_cpu_mask[node]->maskp);
                        free(node_cpu_mask[node]);
                    }

                    free(node_cpu_mask);
                }

                int parse_bitmask(char * line, struct bitmask * mask)
                {
                    char * p = strchr(line, '\n');
                    if (!p)
                        return -1;

                    for (unsigned i(0); p > line; ++i)
                    {
                        char * oldp, * endp;
                        oldp = p;
                        if (*p == ',')
                            --p;

                        while (p > line && *p != ',')
                            --p;

                        if (p > line && sizeof(unsigned long) == 8)
                        {
                            oldp--;
                            memmove(p, p+1, oldp-p+1);
                            while (p > line && *p != ',')
                                --p;
                        }

                        if (*p == ',')
                            p++;

                        if (i >= CPU_LONGS(mask->size))
                            return -1;

                        mask->maskp[i] = strtoul(p, &endp, 16);
                        if (endp != oldp)
                            return -1;

                        p--;
                    }

                    return 0;
                }


                bool create_node_cpu_mask()
                {
                    bool success = true;
                    size_t len(0);
                    char fn[64];
                    FILE * f;

                    for (unsigned node(0) ; node < num_nodes ; ++node)
                    {
                        node_cpu_mask[node] = (bitmask *) malloc(sizeof(struct bitmask *));
                        node_cpu_mask[node]->size = mask_size;
                        node_cpu_mask[node]->maskp = (unsigned long *) calloc(longsperbits(mask_size), sizeof(unsigned long));

                        char * line = NULL;

                        sprintf(fn, "/sys/devices/system/node/node%d/cpumap", node);
                        f = fopen(fn, "r");
                        if (!f || getdelim(&line, &len, '\n', f) < 1)
                            success = false; // SysFS is not available or invalid

                        if (f)
                            fclose(f);

                        if (line && (parse_bitmask(line, node_cpu_mask[node]) < 0))
                            success = false; // Cannot parse bitmask

                        free(line);
                    }

                    return success;
                }

                unsigned node_of_lpu(unsigned lpu)
                {
                    unsigned ret(0);

                    for (unsigned node(0); node < num_nodes; node++)
                    {
                        if (node_cpu_mask[node]->maskp[lpu / bitsperlong] >> (lpu % bitsperlong) & 1)
                        {
                            ret = node;
                            break;
                        }
                    }

                    return ret;
                }
            };

            static unsigned retrieve_num_nodes()
            {
                DIR * d;
                struct dirent * de;

                unsigned num_nodes(0);

                d = opendir("/sys/devices/system/node");
                if (!d)
                    num_nodes = 1;
                else
                {
                    while ((de = readdir(d)) != NULL)
                    {
                        if (0 == strncmp(de->d_name, "node", 4))
                            num_nodes++;
                    }
                    closedir(d);
                }
                return num_nodes;
            }

        static unsigned num_nodes()
        {
            return intern::retrieve_num_nodes();
        }

        static unsigned * cpu_to_node_array(unsigned num_nodes, unsigned num_lpus)
        {
            unsigned * result;

            if (num_nodes < 2)
            {
                return NULL;
            }

            // if num_nodes >= 2 then num_lpus should be >=2 too...
            result = new unsigned[num_lpus];
            intern::NumaInfo info;
            bool success = info.init(num_nodes);

            if (! success)
                return NULL;

            for (unsigned i(0) ; i < num_lpus ; ++i)
            {
                result[i] = info.node_of_lpu(i);
            }

            info.cleanup();

            return result;
        }
        }


    }
}

#endif
