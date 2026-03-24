#!/usr/bin/env bash
# exit immediately on error (-e), treat unset variables as errors (-u),
# and propagate pipe failures (-o pipefail) to catch silent failures in piped commands
set -euo pipefail

# list of all dimensionless cavity distances d/R to simulate
# these correspond to the full set of d/R values used in the paper (Figures 2 and 3)
# spacing is geometrically graded: finer near d/R=1 (close interaction) and coarser at large d/R
yc_values=(1.10 1.25 1.41 1.60 1.82 2.06 2.33 2.65 3.00 3.37 3.78 4.24 4.76 5.35 6.00 6.63 7.33 8.11 8.96 9.91 10.95 12.11 13.39 14.80 16.36 18.09 20.00)

# loop over all d/R values and launch a complete simulation for each
for yc in "${yc_values[@]}"; do
  echo "=== Running y_c=${yc} ==="
  # call run_sim.sh with fixed parameters for all batch runs:
  #   mesh_quality = Finest   → highest accuracy (h_min = 0.002R near cavity)
  #   Delta_Vbarf  = 3374     → large final volume expansion to approach the cavitation limit
  #   steps        = 15000    → fine load increments for numerical stability at large deformations
  #   BC           = fixed    → cavity-interface interaction (rigid wall at bottom)
  # to run cavity-cavity simulations instead, change "fixed" to "symm"
  ./run_sim.sh "$yc" Finest 3374 15000 fixed
done

# -------------------------------------------------------
# Reference: run_sim.sh argument order
# ./run_sim.sh <y_c> <mesh_quality> <Delta_Vbarf> <steps> <BC: symm|fixed>
#
# y_c          : d/R values from the yc_values array above
# mesh_quality : Coarse | Medium | Fine | Finest
# Delta_Vbarf  : (V - V0) / V0 — dimensionless final cavity volume expansion ratio
# steps        : number of incremental load steps in the FEniCS Newton solver
# BC           : "symm" → cavity-cavity interaction (symmetry plane at bottom)
#                "fixed" → cavity-interface interaction (rigid wall at bottom)
# -------------------------------------------------------