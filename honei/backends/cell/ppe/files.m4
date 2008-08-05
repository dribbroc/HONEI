dnl vim: set ft=m4 et :
dnl This file is used by Makefile.am.m4. You should
dnl use the provided autogen.bash script to do all the hard work.
dnl
dnl This file is used to avoid having to make lots of repetitive changes in
dnl Makefile.am every time we add a source or test file. The first parameter is
dnl the base filename with no extension; later parameters can be `hh', `cc',
dnl `test', `impl', `testscript'. Note that there isn't much error checking done
dnl on this file at present...

add(`spe',                            `hh', `cc')
add(`spe_error',                      `hh', `cc')
add(`spe_event',                      `hh', `cc')
add(`spe_instruction',                `hh', `cc')
add(`spe_kernel',                     `hh', `cc', `test')
add(`spe_kernel_manager',             `hh', `cc')
add(`spe_manager',                    `hh', `cc', `test')
add(`spe_transfer_list',              `hh', `cc', `test')

