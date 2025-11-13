#!/bin/bash
set -euo pipefail

###############################################
# Step 0 — Fix PDB residue names for CpHMD
###############################################

RAW_PDB="myc_max_rcsb_rotated.pdb"
PDB_IN="myc_max_rcsb_rotated_cphmd.pdb"
BASE="myc_max_rcsb_rotated_cphmd"
TARGET_M=0.150

echo "=== Step 0: Fixing PDB for CpHMD ==="

if [[ ! -f "$RAW_PDB" ]]; then
    echo "Error: Input file '$RAW_PDB' not found."
    exit 1
fi

sed -E '
    s/\bASP\b/AS2/g;
    s/\bGLU\b/GL2/g;
    s/\bHIS\b/HIP/g;
' "$RAW_PDB" > "$PDB_IN"

echo "Created CpHMD-ready PDB: $PDB_IN"
echo ""

###############################################
# Step 1 — Initial solvation (no ions)
###############################################

echo "=== Step 1: Initial solvation (no ions) ==="

cat > tleap_initial.in <<'TLEAP'
source leaprc.protein.ff19SB
source leaprc.gaff
source leaprc.water.opc

set default PBradii mbondi

loadamberparams frcmod.phmd
loadoff phmd.lib

model = loadpdb MY_PDB
solvateoct model OPCBOX 15 0.75
charge model

savepdb model TMP_SOLV.pdb
saveamberparm model TMP.prmtop TMP.rst7
quit
TLEAP

# Insert actual PDB name
sed -i "s|MY_PDB|${PDB_IN}|g" tleap_initial.in

tleap -f tleap_initial.in > tleap_initial.log

###############################################
# Step 2 — Count waters + net charge
###############################################

echo "=== Step 2: Count waters + charge ==="

# Count waters added by solvation
Nw=$(grep -E "Added [0-9]+ residues" tleap_initial.log | tail -1 | grep -oE '[0-9]+')
if [[ -z "$Nw" ]]; then
    echo "ERROR: Could not parse water count" >&2
    exit 1
fi

# Parse unperturbed charge
k=$(grep -i "unperturbed charge" tleap_initial.log | tail -1 | grep -oP '\([^)]+\)' | tr -d '() ')
if [[ -z "$k" ]]; then
    echo "ERROR: Could not parse unperturbed charge" >&2
    exit 1
fi

echo "Parsed: num_water = $Nw   unperturbed_charge = $k"

###############################################
# Step 3 — Compute ions using voxel volume
###############################################

echo "=== Step 3: Compute ions (target ${TARGET_M} mM) ==="

python compute_solvent_volume_and_ions.py TMP.prmtop TMP.rst7 $TARGET_M

if [[ ! -f tleap_ion_counts.txt ]]; then
    echo "ERROR: Missing tleap_ion_counts.txt" >&2
    exit 1
fi

echo "Ion counts:"
cat tleap_ion_counts.txt
echo ""

###############################################
# Step 4 — Final solvation & ion placement
###############################################

echo "=== Step 4: Final solvation + add ions ==="

cat > tleap_final.in <<TLEAP
source leaprc.protein.ff19SB
source leaprc.gaff
source leaprc.water.opc

set default PBradii mbondi

loadamberparams frcmod.phmd
loadoff phmd.lib

model = loadpdb ${PDB_IN}
solvateoct model OPCBOX 15 0.75

$(cat tleap_ion_counts.txt)

charge model

savepdb model ${BASE}_solvated.pdb
saveamberparm model ${BASE}.prmtop ${BASE}.rst7
quit
TLEAP

tleap -f tleap_final.in > tleap_final.log

###############################################
# Step 5 — Sanity check of actual concentration
###############################################

echo "=== Step 5: SANITY CHECK — compute salt concentration again ==="

python compute_solvent_volume_and_ions.py ${BASE}.prmtop ${BASE}.rst7 | tee sanity_check.log

echo ""
echo "=== Done ==="
echo "Outputs:"
echo "  ${BASE}_solvated.pdb"
echo "  ${BASE}.prmtop"
echo "  ${BASE}.rst7"
echo "Sanity check saved → sanity_check.log"

