dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`general', `assertion',                 `hh', `cc', `test')
add(`general', `condition_variable',        `hh', `cc')
add(`general', `exception',                 `hh', `cc')
add(`hdf5',    `hdf5',                      `hh', `cc', `test')
add(`general', `lock',                      `hh', `cc')
add(`general', `log',                       `hh', `cc')
add(`general', `memory_backend',            `hh', `cc')
add(`cell',    `memory_backend_cell',       `hh', `cc', `test')
add(`gpu',     `memory_backend_gpu',        `hh', `cc')
add(`general', `memory_manager',            `hh', `cc', `test')
add(`general', `mutex',                     `hh', `cc')
add(`general', `pool_task',                 `hh',       `test')
add(`general', `pool_thread',               `hh', `cc', `test')
add(`general', `shared_array',              `hh', `cc', `test')
add(`cell',    `spe',                       `hh', `cc')
add(`cell',    `spe_event',                 `hh', `cc')
add(`cell',    `spe_instruction',           `hh', `cc')
add(`cell',    `spe_kernel',                `hh', `cc', `test')
add(`cell',    `spe_manager',               `hh', `cc', `test')
add(`general', `stringify',                 `hh')
add(`general', `sync_point',                `hh')
add(`general', `tags',                      `hh', `cc')
add(`general', `thread',                    `hh', `cc', `test')
add(`general', `thread_pool',               `hh',       `test')
add(`general', `transfer',                  `hh', `cc')
add(`general', `type_traits',               `hh', `cc')
add(`general', `worker',                    `hh', `cc', `test')
add(`general', `wrapper',                   `hh')

