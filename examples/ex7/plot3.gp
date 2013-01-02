
#set term postscript enhanced color
set term post eps
set output 'out1.eps'

set datafile separator ' '

set format y "%.2e"

set multiplot layout 2,1;                          # get into multiplot mode

set ylabel 'Stress [MPa]' 
set xlabel 'Strain' 
unset key
plot 'constitutive_concrete.txt' u 1:2 w l


unset key
plot 'constitutive_steel.txt' u 1:2 w l


unset multiplot                         # exit multiplot mode

