# need to construct mini_?.min files
pmemd.cuda -O -i heating.mdin -p myc_max_rcsb_rotated_cphmd.prmtop -c myc_max_rcsb_rotated_cphmd_min1.rst7 -r myc_max_rcsb_rotated_cphmd_heating.rst7 -ref myc_max_rcsb_rotated_cphmd_min1.rst7 -phmdparm ff19sb_pme.parm -phmdin prod.phmdin -o heating.mdout -x myc_max_rcsb_rotated_cphmd_heating.nc
