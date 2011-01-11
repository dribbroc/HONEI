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
add(`general', `barrier',                        `hh', `cc')
add(`general', `benchmark_info',                 `hh')
add(`general', `condition_variable',             `hh', `cc')
add(`general', `configuration',                  `hh', `cc', `test')
add(`general', `exception',                      `hh', `cc')
add(`general', `file_to_string',                 `hh')
add(`hdf5',    `hdf5',                           `hh', `cc', `test')
add(`general', `instantiation_policy',           `hh', `impl')
add(`hdf5',    `kpnetcdffile',                   `hh', `cc', `test')
add(`hdf5',    `kpnetcdf_types',                 `hh')
add(`general', `lock',                           `hh', `cc')
add(`general', `log',                            `hh', `cc')
add(`general', `memory_arbiter',                 `hh', `impl', `cc', `test')
add(`general', `memory_backend_base',            `hh')
add(`general', `memory_backend',                 `hh', `cc', `test')
add(`general', `memory_pool',                    `hh', `cc', `test')
add(`general', `mutex',                          `hh', `cc')
add(`hdf5',    `netcdf_datatypes',               `hh')
add(`general', `operation_wrapper',              `hh')
add(`general', `partitioner',                    `hh', `cc', `test')
add(`general', `private_implementation_pattern', `hh', `impl')
add(`general', `profiler',                       `hh', `cc', `test')
add(`general', `shared_array',                   `hh', `impl', `cc', `test')
add(`general', `stringify',                      `hh')
add(`general', `string_tokenizer',               `hh')
add(`general', `sync_point',                     `hh')
add(`general', `tags',                           `hh', `cc')
add(`general', `ticket',                         `hh', `cc')
add(`general', `time_stamp',                     `hh', `test')
add(`general', `thread',                         `hh', `cc', `test')
add(`general', `type_traits',                    `hh', `cc', `test')
add(`general', `tr1_boost',                      `hh')
add(`general', `unittest',                       `hh', `cc', `test')
add(`general', `wrapped_forward_iterator',       `hh', `fwd', `impl')
