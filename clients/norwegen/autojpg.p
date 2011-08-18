set output "FILENAME.jpg"
set palette rgbformulae 23,28,3
set terminal jpeg
set notitle
set view 60,30
unset xtics
unset ytics
unset ztics
set nokey
unset colorbox
set border 4095

set xrange[*:*]
set yrange[*:*]
set zrange[*:*]

splot\
"FILENAME" notitle with pm3d
