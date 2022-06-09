set size ratio 0.8

aqua = "#00FFFF"; azure = "#F0FFFF"
aliceblue = "#F0F8FF"

set key font ",18"

set grid ytics mytics  # draw lines for each ytics and mytics
set grid xtics mytics
set mytics 1           # set the spacing for the mytics
set grid               # enable the grid

set style fill solid 1.00 border 0

set xlabel "t/T" font ",18"

set ylabel "<Logical_{Z}>" font ",18"
# set ylabel offset -2,0

set xtics font ",10"
set ytics font ",10"

set xrange[0:50]
set yrange[-1:1]





set title "" font ",24"



plot 'nodecoder_N=12_MaxDim=100_ns=200.dat' using 1:2 with lines ls 1 lw 2 t "200",\
'nodecoder_N=12_MaxDim=100_ns=200.dat' using 1:3 with lines ls 1 lw 2 t "",\
'nodecoder_N=12_MaxDim=100_ns=400.dat' using 1:2 with lines ls 3 lw 2 t "400",\
'nodecoder_N=12_MaxDim=100_ns=400.dat' using 1:3 with lines ls 3 lw 2 t "",\
'nodecoder_N=12_MaxDim=500.dat' using 1:2 with lines ls 4 lw 2 t "500",\
'nodecoder_N=12_MaxDim=500.dat' using 1:3 with lines ls 4 lw 2 t "",\
'nodecoder_N=12_MaxDim=200.dat' using 1:2 with lines ls 5 lw 2 t "200",\
'nodecoder_N=12_MaxDim=200.dat' using 1:3 with lines ls 5 lw 2 t "",\
'nodecoder_dx=3.dat' using 1:2 with lines ls 2 lw 2 t "normal",\
'nodecoder_dx=3.dat' using 1:3 with lines ls 2 lw 2 t "",\


# plot 'nodecoder_N=20_MaxDim=200.dat' using 1:2 with lines ls 1 lw 2 t "N sample = 4000",\
# 'nodecoder_N=20_MaxDim=200.dat' using 1:3 with lines ls 1 lw 2 t "",\
# 'nodecoder_dx=5.dat' using 1:2 with lines ls 2 lw 2 t "N sample = 8000",\
# 'nodecoder_dx=5.dat' using 1:3 with lines ls 2 lw 2 t "",\

