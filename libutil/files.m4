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
add(`general', `lock',                      `hh', `cc')
add(`general', `log',                       `hh', `cc')
add(`general', `memory_backend',            `hh')
add(`gpu',     `memory_backend_gpu',        `hh', `cc')
add(`general', `memory_manager',            `hh', `cc', `test')
add(`general', `mutex',                     `hh', `cc')
add(`general', `shared_array',              `hh', `cc')
add(`cell',    `spe_manager',               `hh', `cc', `test')
add(`general', `stringify',                 `hh')
add(`general', `tags',                      `hh', `cc')
