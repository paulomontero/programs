set terminal pdf enhanced color font ',12'
set output "pbh_h12_rho.pdf"

set border linewidth 2

set xlabel '{/Symbol r}		[g {cm}^{-3}]'
set ylabel 'H_{12}'

set xrange [1e6:1e10]
#set yrange [0:1.0]
#set xlabel 'k 		[{Mpc}^{-1}]'
#set xlabel 'z'
#set ylabel '{/Symbol D}_{m,{/Symbol y}}^{3D} / {/Symbol D}_{m}^{3D}'
#set ylabel 'P_{m,{/Symbol y}}^{1D} / P_{m}^{1D}'
#set ylabel '{/Symbol D}_{m,x_{HI}}'
#set ylabel '{/Symbol D}_{m,{/Symbol y}}'
set logscale x
#set logscale y



set style line 2 lt 1 lc rgb "dark-green" lw 4 pt 6 ps 1 #--green


set key center left

#set grid
set tics scale 1.5

#set arrow from 0.215031, graph(0,0) to 0.215031, graph(1,1) linetype 0 nohead 

plot 'plot_h12_rho.dat' using 1:2 title 'T_{9} = 7' with lines lt 7 lc rgb "dark-green" lw 3, 'plot_h12_rho.dat' using 1:3 title 'T_{9} = 5' with lines lt 7 lc rgb "navy" lw 3, 'plot_h12_rho.dat' using 1:4 title 'T_{9} = 3' with lines lt 7 lc rgb "orange" lw 3, 'plot_h12_rho.dat' using 1:5 title 'T_{9} = 1' with lines lt 7 lc rgb "purple" lw 3

#plot 'plot_temp_t.dat' using 3:10 notitle with lp pt 8 ps 0.3 lt 7 lc 7 lw 1
#plot 'plot_temp_max_test.dat' using 3:2 title 'T_{max}' with lines lt 7 lc rgb "dark-green" lw 3, 'plot_temp_final_test.dat' using 3:2 title 'T_{final}' with lines lt 7 lc rgb 'purple' lw 3
#plot 'vz_final.dat' using 2:3 title 'timestep = 1' with lines lt 7 lc rgb "magenta" lw 3

#plot 'new_plot_vz_tstep5001.dat' using 4:3 notitle with lines lt 7 lc rgb "orange" lw 2

#plot 'vz_tstep500.dat' using 4:3 title 'timestep = 500' with lines lt 7 lc rgb "navy" lw 3, 'vz_tstep400.dat' using 4:3 title 'timestep = 400' with lines lt 7 lc rgb "olive" lw 3, 'vz_tstep250.dat' using 4:3 title 'timestep = 250' with lines lt 7 lc rgb "orange" lw 3, 'vz_tstep100.dat' using 4:3 title 'timestep = 100' with lines lt 7 lc rgb "magenta" lw 3
#plot 'vz_tstep500.dat' using 4:3 title 'timestep = 500' with lines lt 7 lc rgb "navy" lw 3


