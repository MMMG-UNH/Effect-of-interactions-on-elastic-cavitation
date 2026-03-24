#!/usr/bin/env bash
# exit immediately on error (-e), treat unset variables as errors (-u),
# and propagate pipe failures (-o pipefail) to catch silent failures in piped commands
set -euo pipefail

# print usage instructions if the wrong number of arguments is provided
usage() {
  echo "Usage: $0 <y_c> <mesh_quality> <Delta_Vbarf> <steps> <BC: symm|fixed>"
  exit 1
}

# enforce that exactly 5 arguments are passed before proceeding
if [[ $# -ne 5 ]]; then
  usage
fi

# y_c: dimensionless cavity center position d/R
#   - for cavity-interface problem: distance from cavity center to rigid interface / cavity radius
#   - for cavity-cavity problem:    half the center-to-center distance / cavity radius
yc=$1              # e.g. 1.10

# mesh_quality: controls element size near the cavity and in the far field
# options: Coarse | Medium | Fine | Finest  (defined in geo_file_creator_2.py)
mq=$2              # e.g. Finest

# Delta_Vbarf: target dimensionless final cavity volume expansion = (V_final - V0) / V0
delta_Vbarf=$3     # e.g. 3374

# steps: number of incremental pseudo-time load steps in the FEniCS solver
steps=$4           # e.g. 15000

# bc: selects the bottom boundary condition type, which determines which problem is solved
#   symm  → symmetry BC (u2=0 only)  → models cavity-cavity interaction (two identical cavities)
#   fixed → fixed BC  (u1=u2=0)      → models cavity-interface interaction (cavity near rigid wall)
bc=$5              # "symm" or "fixed"

# construct the .geo and mesh filenames from the input parameters
meshgeo="2Cavity_yc${yc}_${mq}.geo"
meshname="2Cavity_yc${yc}_${mq}"

# step 1: generate the Gmsh .geo geometry file for this d/R value and mesh quality
echo "→ Generating .geo file..."
python3 geo_file_creator_2.py "$yc" "$mq"

# step 2: move the generated .geo file into the meshes/ subdirectory
echo "→ Moving .geo → meshes/"
mv "$meshgeo" meshes/

# step 3: run Gmsh in 2D batch mode to generate the .msh mesh file from the .geo file
# -2 flag: generate a 2D mesh; -format msh2: use legacy MSH2 format required by meshio
echo "→ Running Gmsh"
cd meshes
gmsh -2 -format msh2 "$meshgeo"
cd ..

# step 4: convert the Gmsh .msh file to XDMF format for FEniCS
# produces two files: <meshname>.xdmf (domain) and facet_<meshname>.xdmf (boundary tags)
echo "→ Converting to XDMF"
cd meshes
python3 ../2d_gmsh_convert.py "$meshname"
cd ..

# step 5: launch the FEniCS boundary value problem solver
# all runs use remote displacement loading (arg 1 = 1) with free horizontal outer boundary (arg 5 = 0)
echo "→ Launching FEniCS with:"
echo "     remote_disp=1"
echo "     Delta_Vbarf=${delta_Vbarf}"
echo "     steps=${steps}"
echo "     save_xdmf=1"
echo "     remote_disp_uxfix=0"
echo "     y_c=${yc}"
echo "     mesh_quality=${mq}"

# dispatch to the correct BVP script based on the boundary condition type
if [[ "$bc" == "symm" ]]; then
  # symmetry BC: cavity-cavity interaction problem (mirror plane at bottom)
  python3 BVP_Bulk_2Cavity_symmBC_2.py 1 "$delta_Vbarf" "$steps" 1 0 "$yc" "$mq"
elif [[ "$bc" == "fixed" ]]; then
  # fixed BC: cavity-interface interaction problem (rigid wall at bottom)
  python3 BVP_Bulk_2Cavity_fixedBC_2.py 1 "$delta_Vbarf" "$steps" 1 0 "$yc" "$mq"
fi

# confirm successful completion with a summary of the run parameters
echo "✅ Completed simulation for: y_c=${yc}, mesh_quality=${mq}, Delta_Vbarf=${delta_Vbarf}, steps=${steps}, BC=${bc}"