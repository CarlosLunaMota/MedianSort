set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 1200, 600 
set xlabel 'Size of the array'
set ylabel 'Time (% w.r.t qsort)'
set xtics  out nomirror
set ytics  out nomirror
set logscale x
set key autotitle columnhead outside

min(x,y) = (x < y) ? x : y
max(x,y) = (x > y) ? x : y

set output 'plot1.png'
set title "{/=20 MedianSort vs QuickSort vs HeapSort}\n(average performance for 10 uniformly random instances of each size)"
set yrange [45:110]

plot 'log.txt' using 1:2  with linespoints lw 3 pt  7 lc rgb "#000000", \
     'log.txt' using 1:3  with linespoints lw 3 pt  7 lc rgb "#119911", \
     'log.txt' using 1:4  with linespoints lw 2 pt  5 lc rgb "#FFAAAA", \
     'log.txt' using 1:5  with linespoints lw 2 pt  7 lc rgb "#FF9999", \
     'log.txt' using 1:6  with linespoints lw 2 pt  9 lc rgb "#FF8888", \
     'log.txt' using 1:7  with linespoints lw 2 pt 11 lc rgb "#FF7777", \
     'log.txt' using 1:8  with linespoints lw 2 pt 13 lc rgb "#FF6666", \
     'log.txt' using 1:9  with linespoints lw 2 pt  5 lc rgb "#FF5555", \
     'log.txt' using 1:10 with linespoints lw 2 pt  7 lc rgb "#FF4444", \
     'log.txt' using 1:11 with linespoints lw 2 pt  9 lc rgb "#FF3333", \
     'log.txt' using 1:12 with linespoints lw 2 pt 11 lc rgb "#FF2222", \
     'log.txt' using 1:13 with linespoints lw 2 pt 13 lc rgb "#FF1111", \
     'log.txt' using 1:14 with linespoints lw 2 pt  5 lc rgb "#AAAAFF", \
     'log.txt' using 1:15 with linespoints lw 2 pt  7 lc rgb "#9999FF", \
     'log.txt' using 1:16 with linespoints lw 2 pt  9 lc rgb "#8888FF", \
     'log.txt' using 1:17 with linespoints lw 2 pt 11 lc rgb "#7777FF", \
     'log.txt' using 1:18 with linespoints lw 2 pt 13 lc rgb "#6666FF", \
     'log.txt' using 1:19 with linespoints lw 2 pt  5 lc rgb "#5555FF", \
     'log.txt' using 1:20 with linespoints lw 2 pt  7 lc rgb "#4444FF", \
     'log.txt' using 1:21 with linespoints lw 2 pt  9 lc rgb "#3333FF", \
     'log.txt' using 1:22 with linespoints lw 2 pt 11 lc rgb "#2222FF"


set output 'plot2.png'
set title "{/=20 MedianSort vs QuickSort vs HeapSort}\n(average performance for 10 uniformly random instances of each size)"

plot 'log.txt' using 1:(min($4,min($5,min($6,min($7,min($8,min($9,min($10,min($11,min($12,$13)))))))))):(max($4,max($5,max($6,max($7,max($8,max($9,max($10,max($11,max($12,$13)))))))))) with filledcurves lc rgb "#CC1111" fs transparent solid 0.15 notitle, \
     'log.txt' using 1:(min($14,min($15,min($16,min($17,min($18,min($19,min($20,min($21,$22))))))))):(max($14,max($15,max($16,max($17,max($18,max($19,max($20,max($21,$21))))))))) with filledcurves lc rgb "#1111CC" fs transparent solid 0.15 notitle, \
     'log.txt' using 1:2  with linespoints lw 3 pt 7 lc rgb "#000000", \
     'log.txt' using 1:3  with linespoints lw 3 pt 7 lc rgb "#119911", \
     'log.txt' using 1:4  with linespoints lw 2 pt 0 lc rgb "#CC1111", \
     'log.txt' using 1:10 with linespoints lw 3 pt 7 lc rgb "#CC1111", \
     'log.txt' using 1:14 with linespoints lw 2 pt 0 lc rgb "#1111CC", \
     'log.txt' using 1:20 with linespoints lw 3 pt 7 lc rgb "#1111CC"
