# need to construct mini_?.min files
pmemd.cuda -O -i heating.mdin -p 1w4h_cphmd.prmtop -c 1w4h_min1.rst7 -r 1w4h_heating.rst7 -ref 1w4h_min1.rst7 -phmdparm ff19sb_pme.parm -phmdin prod.phmdin -o heating.mdout -x 1w4h_heating.nc
