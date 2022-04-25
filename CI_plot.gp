#!/usr/local/Cellar/gnuplot/5.0.5_2/bin/gnuplot
#
#    
#    	G N U P L O T
#    	Version 5.0 patchlevel 5    last modified 2016-10-02
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2016
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal aqua 0 title "Figure 0" size 846,594 font "Times-Roman,14" enhanced solid
# set output

GNUTERM = "aqua"

#
set macros
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # BLUE
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5  # RED
set style line 3 lc rgb 'green' lt 1 lw 2 pt 5       #GREEN


set style line 4 lc rgb 'web-blue' lt 1 lw 2 pt 5      #WEB BLUE
set style line 5 lc rgb 'orange' lt 1 lw 2 pt 5       #ORANGE
set style line 6 lc rgb 'dark-green' lt 1 lw 2 pt 5       #DARK GREEN

set style line 7 lc rgb 'violet' lt 1 lw 2 pt 5      #VIOLET
#
BLUE = '1'
RED = '2'
GREEN = '3'
WEBBLUE = '4'
ORANGE = '5'
DARKGREEN = '6'
VIOLET = '7'
#
set grid
#
set term aqua
set key
set title "CI percolation " font "sans,16"
set xlabel "fraction CI nodes removed" font "Helvetica, 14"
set ylabel "GC" font "Helvetica, 14"
set yrange[0:1.01]
set xrange[0:0.7]
set xtics font "Helvetica, 14" 0,0.05
set mxtics 5
set ytics font "Helvetica, 14" 0,0.05
set mytics 5
#
plot 'CI_NoN_L_0.txt' u 2:3 title 'CI 0' w lp ls @BLUE
replot 'CI_NoN_L_1.txt' u 2:3 title 'CI 1' w lp ls @RED
replot 'CI_NoN_L_2.txt' u 2:3 title 'CI 2' w lp ls @GREEN
replot 'CI_NoN_L_3.txt' u 2:3 title 'CI 3' w lp ls @ORANGE
#replot 'CI_L_4.txt' u 1:3 title 'CI 4' w lp ls @DARKGREEN


#    EOF

