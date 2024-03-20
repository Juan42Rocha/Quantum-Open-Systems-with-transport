set term wxt persist


set grid
set k l
set log y


pl "nu0.0001-1.dat" u 1:2 w lp
repl "nu0.001-1.dat" u 1:2 w lp
repl "nu0.1-2.dat" u 1:2 w lp
repl "nu0.2-2.dat" u 1:2 w lp
repl "nu0.5-2.dat" u 1:2 w lp
repl "nu0.7-2.dat" u 1:2 w lp
