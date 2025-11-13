# need to construct mini_?.min files
pmemd.cuda -O -i equil_1.1.mdin -p 1w4h_cphmd.prmtop -c 1w4h_heating.rst7 -r 1w4h_equil_1.1.rst7 -ref 1w4h_heating.rst7 -phmdparm ff19sb_pme.parm -phmdin prod.phmdin -phmdstrt phmdrestrt -o equil_1.1.mdout -x 1w4h_equil1.nc
