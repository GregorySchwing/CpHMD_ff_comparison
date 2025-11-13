# need to construct mini_?.min files
pmemd.cuda -O -i minimization1.mdin -p 1w4h_cphmd.prmtop -c 1w4h_cphmd.rst7 -r 1w4h_min1.rst7 -o 1w4h_min_1.out -ref  1w4h_cphmd.rst7 -x 1w4h_min1.nc
