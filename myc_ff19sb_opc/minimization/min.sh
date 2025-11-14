# need to construct mini_?.min files
cp ../myc_max_rcsb_rotated_cphmd.rst7 .
cp ../myc_max_rcsb_rotated_cphmd.prmtop .
pmemd.cuda -O -i minimization1.mdin -p myc_max_rcsb_rotated_cphmd.prmtop -c myc_max_rcsb_rotated_cphmd.rst7 -r myc_max_rcsb_rotated_cphmd_min1.rst7 -o myc_max_rcsb_rotated_cphmd_min1.out -ref  myc_max_rcsb_rotated_cphmd.rst7 -x myc_max_rcsb_rotated_cphmd_min1.nc
