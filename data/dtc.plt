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





# set title "1 sample" font ",24"
#
# plot 'data_ns=1_nt=1000.dat' using 1:2 with lines lc 1 lw 2 t "Without decoder",\
# 'data_ns=1_nt=1000.dat' using 1:3 with lines lc 2 lw 2 t "Without decoder",\
# 'data_decoder_thd=0.100000_ns=1_nt=1000.dat' using 1:2 with lines lc 3 lw 2 t "With decoder",\
# 'data_decoder_thd=0.100000_ns=1_nt=1000.dat' using 1:3 with lines lc 4 lw 2 t "With decoder",\


set title "" font ",24"

#
# plot 'nodecoder_dx=3_ns=50_nt=100.dat' using 1:2 with lines ls 1 lw 2 t "1",\
# 'nodecoder_dx=3_ns=50_nt=100.dat' using 1:3 with lines ls 1 lw 2 t "",\
# 'decoder_dx=3_ns=50_nt=50.dat' using 1:2 with lines ls 2 lw 2 t "2",\
# 'decoder_dx=3_ns=50_nt=50.dat' using 1:3 with lines ls 2 lw 2 t "",\


plot 'nodecoder_N=12_MaxDim=100.dat' using 1:2 with lines ls 1 lw 2 t "maxdim = 100",\
'nodecoder_N=12_MaxDim=100.dat' using 1:3 with lines ls 1 lw 2 t "",\
'nodecoder_N=12.dat' using 1:2 with lines ls 2 lw 2 t "maxdim = 1000",\
'nodecoder_N=12.dat' using 1:3 with lines ls 2 lw 2 t "",\
'nodecoder_dx=3.dat' using 1:2 with lines ls 3 lw 2 t "normal",\
'nodecoder_dx=3.dat' using 1:3 with lines ls 3 lw 2 t "",\


