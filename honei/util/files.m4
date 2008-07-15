dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`general', `assertion',                      `hh', `cc', `test')
add(`general', `attributes',                     `hh')
add(`general', `benchmark_info',                 `hh')
add(`general', `condition_variable',             `hh', `cc')
add(`general', `configuration',                  `hh', `cc', `test')
add(`general', `exception',                      `hh', `cc')
add(`hdf5',    `hdf5',                           `hh', `cc', `test')
add(`general', `instantiation_policy',           `hh', `impl')
add(`general', `lock',                           `hh', `cc')
add(`general', `log',                            `hh', `cc')
add(`general', `memory_arbiter',                 `hh', `impl', `cc')
add(`general', `memory_backend_base',            `hh')
add(`general', `memory_backend',                 `hh', `cc')
add(`general', `mutex',                          `hh', `cc')
add(`general', `partitioner',                    `hh', `cc', `test')
add(`general', `pool_task',                      `hh',       `test')
add(`general', `pool_thread',                    `hh', `cc', `test')
add(`general', `private_implementation_pattern', `hh', `impl')
add(`general', `profiler',                       `hh', `cc', `test')
add(`general', `shared_array',                   `hh', `impl', `cc', `test')
add(`general', `stringify',                      `hh')
add(`general', `sync_point',                     `hh')
add(`general', `tags',                           `hh', `cc')
add(`general', `time_stamp',                     `hh', `test')
add(`general', `thread',                         `hh', `cc', `test')
add(`general', `thread_pool',                    `hh',       `test')
add(`general', `type_traits',                    `hh', `cc', `test')
add(`general', `worker',                         `hh', `cc', `test')
add(`general', `wrapper',                        `hh')

