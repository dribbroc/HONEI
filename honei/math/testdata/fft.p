set output "fft.eps"
set terminal postscript eps enhanced color
#set size 0.75,0.75
set notitle
#set yrange [0.24:0.25]
#set xlabel "x"
#set ylabel "h + b"
#set xtics ("0" 0, "1" 100, "2" 199)
#set line width 2
plot\
"fft.dat" with points lt 2 lw 2

