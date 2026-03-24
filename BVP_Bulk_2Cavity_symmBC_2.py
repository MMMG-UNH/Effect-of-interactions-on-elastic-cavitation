"""
BVP of axisymmetric spherical cavity expansion in bulk of incompressible material with surface tension
(set to zero here) at cavity surface.
Cavity volume can be controlled by applied displacement on outer boundary
or cavity surface.

This script solves the cavity-cavity interaction problem: two identical spherical cavities
in an infinite medium. The symmetry plane between the two cavities is modeled as the bottom
boundary with a symmetry BC (u2=0 only — free to slide tangentially), which is equivalent
to placing a mirror image cavity below. This is in contrast to the fixedBC version, which
enforces both u1=0 and u2=0 at the bottom to simulate a rigid interface.

How to run:
    python3 BVP_Bulk_2Cavity_symmBC_2.py A B C D E F G
    A: remote_disp     — 1 = apply displacement at outer (top) boundary, 0 = apply at cavity surface
    B: Delta_Vbarf     — target dimensionless volume expansion (V-V0)/V0
    C: steps           — number of incremental load steps
    D: save_xdmf       — 1 = save displacement/stress fields to .xdmf for ParaView
    E: remote_disp_uxfix — 1 = also fix horizontal displacement at outer boundary, 0 = leave free
    F: y_c             — dimensionless cavity half-spacing d/R (half the center-to-center distance)
    G: mesh_quality    — Coarse / Medium / Fine / Finest

Example:
    python3 BVP_Bulk_2Cavity_symmBC_2.py 1 3374 15000 1 0 1.10 Finest

Output (.mat file columns):
    [0] remote displacement applied at top boundary
    [1] current cavity volume V
    [2] bulk strain energy
    [3] deformed inner (cavity) surface area
    [4] remote pressure (average normal traction on outer boundary)
    [5] inner pressure (average normal traction on cavity surface)
"""

# FEniCS core library
from dolfin import *
import numpy as np
import sys
from datetime import datetime

# set to 1 to save results in .pickle format (in addition to or instead of .mat)
store_data_pickle = 0
# set to 1 to save results in .mat format (readable by MATLAB and scipy)
store_data_mat = 1

# set to 1 to write displacement and stress fields to XDMF for ParaView visualization
save_xdmf = int(sys.argv[4])

# set to 1 to print cavity volume at each step and compare with analytical estimate
print_volume_check = 1
# set to 1 to print deformed coordinates at key points for detailed verification
detailed_check = 1

if store_data_pickle == 1:
    import pickle
if store_data_mat == 1:
    import scipy.io

# import the constitutive model class and result-writing utility
from BulkandSurfaceElasticity_edited import *

# increase quadrature accuracy and enable C++ optimization for form compilation
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["form_compiler"]["cpp_optimize"] = True

# 1 = displacement applied at outer boundary (remote loading), 0 = applied at cavity surface
remote_displacement = int(sys.argv[1])
# y_c: half the center-to-center distance between the two cavities, normalized by cavity radius R
# i.e. y_c = d/R where d is the distance from the symmetry plane to the cavity center
y_c = float(sys.argv[6])
# mesh quality label — must match a key used in geo_file_creator_2.py
quality = sys.argv[7].lower()

# simulation type 2 = axisymmetric (the r-z plane is the simulation domain)
simulation_type = 2

# use incompressible neo-Hookean model
elastic_incompressible = 1
bulk_elastic_model = 1
surface_elastic_model = 1

# shear modulus (set to 1 so all stresses are dimensionless, normalized by mu)
Gshear = 1
# large penalty bulk modulus to enforce near-incompressibility (J ≈ 1)
Kbulk = Gshear*1E12
bulkparams = [Gshear, Kbulk]

# surface energy coefficient at cavity surface (set to 0 = no surface tension)
gambar = 0
surfaceparams = [gambar]

# construct mesh filename from y_c and quality (must match output of geo_file_creator_2.py)
meshname_string = "2Cavity_yc{0:.2f}_{1}".format(y_c, quality)

# instantiate the constitutive model object
BVP = Bulk_SurfaceElasticity(elastic_incompressible, bulk_elastic_model,
                              surface_elastic_model, bulkparams, surfaceparams)

# physical group tags as defined in the .geo file (must match geo_file_creator_2.py)
# Physical Curve("left-top", 10)
# Physical Curve("cavity", 11)
# Physical Curve("left-bottom", 12)
# Physical Curve("bottom", 13)
# Physical Curve("right", 14)
# Physical Curve("top", 15)
Outer_boundary_surface = 15     # top edge — where remote displacement is applied
Inner_boundary_surface = 11     # cavity arc — where surface energy is applied
XZero_boundary_surface = 10     # left edge above cavity — symmetry axis (u1=0)
XZero_boundary_surface2 = 12    # left edge below cavity — symmetry axis (u1=0)
XZero_boundary_surface3 = 14    # right edge — zero horizontal displacement (u1=0)
# bottom edge — symmetry plane between the two cavities:
# only u2=0 is enforced (no normal penetration), u1 is left FREE to allow tangential slip
# this correctly models the symmetry plane of the two-cavity problem
YZero_boundary_surface = 13

# apply surface energy on the cavity surface only
surface_energy_boundaries = Inner_boundary_surface
BVP.setup(meshname_string, simulation_type, surface_energy_boundaries)

# domain size parameters: BbyA = B/A = outer radius / cavity radius = 100
BbyA = 100
# A_param: undeformed cavity radius
A_param = 0.5
# B_param: outer domain radius (= 100 * 0.5 = 50, consistent with L=50 in geo file)
B_param = BbyA * A_param

# target dimensionless volume expansion ratio: (V_final - V0) / V0
Delta_Vbarf = float(sys.argv[2])

# pseudo-time stepping: t goes from 0 to T_f in 'steps' equal increments
t    = 0.0
dt   = 1.0
steps = int(sys.argv[3])
T_f = steps * dt

# flag to also fix horizontal (u1) displacement at the outer boundary
remote_disp_uxfix = int(sys.argv[5])

# total remote vertical displacement at outer boundary corresponding to Delta_Vbarf
# derived from volume conservation for a spherical cavity in a large domain
uy_remote = (4./3.) * B_param * (Delta_Vbarf / BbyA**3.0)
# deformed cavity radius at final step (from volume conservation: V/V0 = (a/A)^3)
a_tot = A_param * pow((1 + Delta_Vbarf), 1./3.)

# set total displacement magnitude and which boundary it is applied on
if(remote_displacement == 1):
    disp_tot = uy_remote
    applied_disp_surface = Outer_boundary_surface
else:
    # if loading from cavity surface, total radial expansion = a_tot - A
    disp_tot = a_tot - A_param
    applied_disp_surface = Inner_boundary_surface

# define time-dependent displacement boundary conditions
# for cavity surface loading: radial displacement in (r, z) components
if(remote_displacement == 0):
    dispx = Expression(("disp_tot*(t/T_f)*x[0]/sqrt(x[0]*x[0] + x[1]*x[1])"),
                        disp_tot=disp_tot, t=0.0, T_f=T_f, degree=2)
    dispy = Expression(("disp_tot*(t/T_f)*x[1]/sqrt(x[0]*x[0] + x[1]*x[1])"),
                        disp_tot=disp_tot, t=0.0, T_f=T_f, degree=2)
else:
    # for remote loading: uniform vertical displacement at the top boundary
    dispx = 0.0
    dispy = Expression(("disp_tot*(t/T_f)"), disp_tot=disp_tot, t=0.0, T_f=T_f, degree=2)

# --- Dirichlet boundary conditions ---
# u2=0 at bottom: zero vertical displacement (symmetry plane — no normal penetration)
bcs_1 = DirichletBC(BVP.ME.sub(0).sub(1), 0, BVP.facets, YZero_boundary_surface)
# u1=0 on left edge above cavity: symmetry axis condition
bcs_2 = DirichletBC(BVP.ME.sub(0).sub(0), 0, BVP.facets, XZero_boundary_surface)
# u1=0 on left edge below cavity: symmetry axis condition (split into two tags in the .geo file)
bcs_6 = DirichletBC(BVP.ME.sub(0).sub(0), 0, BVP.facets, XZero_boundary_surface2)
# applied vertical displacement at either the top or cavity boundary (ramps from 0 to disp_tot)
bcs_3 = DirichletBC(BVP.ME.sub(0).sub(1), dispy, BVP.facets, applied_disp_surface)
# applied horizontal displacement (only used for cavity surface loading)
bcs_4 = DirichletBC(BVP.ME.sub(0).sub(0), dispx, BVP.facets, applied_disp_surface)
# u1=0 on right edge: prevents rigid body horizontal motion
bcs_5 = DirichletBC(BVP.ME.sub(0).sub(0), 0, BVP.facets, XZero_boundary_surface3)
# NOTE: unlike fixedBC, there is NO u1=0 constraint at the bottom boundary —
# the symmetry plane is free to slide tangentially, correctly representing
# the mid-plane between two identical cavities expanding symmetrically

# assemble the active set of boundary conditions based on loading mode
if(remote_displacement == 1 and remote_disp_uxfix == 0):
    # remote loading with free horizontal motion at top: fix axis, bottom (u2 only), right
    bcs = [bcs_1, bcs_2, bcs_3, bcs_5, bcs_6]
else:
    bcs = [bcs_1, bcs_2, bcs_3, bcs_4, bcs_5]

# =========================================================
# Set up the nonlinear variational problem and solver
# =========================================================

# define the nonlinear problem: residual L, solution w, BCs, and Jacobian a
stressProblem = NonlinearVariationalProblem(BVP.L, BVP.w, bcs, J=BVP.a)

# Newton-Raphson nonlinear solver
solver = NonlinearVariationalSolver(stressProblem)

# solver parameters
prm = solver.parameters
prm['nonlinear_solver'] = 'newton'
# MUMPS: a direct parallel sparse solver — robust for large nonlinear elasticity problems
prm['newton_solver']['linear_solver'] = "mumps"
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['relative_tolerance'] = 1e-6
# maximum Newton iterations per load step before declaring non-convergence
prm['newton_solver']['maximum_iterations'] = 30

# suppress FEniCS solver output (log level 30 = WARNING only)
set_log_level(30)

# =========================================================
# Start incremental loading analysis
# =========================================================
print("------------------------------------")
print("Simulation Start")
print("------------------------------------")

print("Simulation Parameters:")
print(f"Simulation Type: {simulation_type}")
print(f"Elastic Incompressible: {elastic_incompressible}")
print(f"Bulk Elastic Model: {bulk_elastic_model}")
print(f"Surface Elastic Model: {surface_elastic_model}")
print(f"Bulk Elasticity Parameters: Gshear = {Gshear}, Kbulk = {Kbulk}")
print(f"Surface Elasticity Parameter: gambar = {gambar}")
print(f"Mesh Name: {meshname_string}")
print(f"Final Volume Expansion Ratio: Delta_Vbarf = {Delta_Vbarf}")
print(f"Number of Steps: {steps}")
print(f"Remote Displacement: {remote_displacement}")
print(f"Remote Displacement Ux Fixed: {remote_disp_uxfix}")

startTime = datetime.now()

# base folder name for saving results
folder_name = "SphericalAxisym2Cavity_SqOuter"

if(remote_displacement == 1):
    folder_name += "RemoteDisp"
    if(remote_disp_uxfix == 0):
        folder_name += "Uxunfixed"

# output file path: encodes key simulation parameters in the filename for traceability
# "_symmBC" suffix distinguishes this output from the fixedBC version
savename = ("results/{0}/DeltaVfbar{1:.2f}_Steps{2:d}_gambar{3:.2f}"
            "_dR_{4:.2f}_symmBC"
            .format(folder_name, Delta_Vbarf, int(steps), gambar, y_c))

if(save_xdmf == 1):
    savenamexdmf = savename + ".xdmf"
    file_results = XDMFFile(savenamexdmf)
    # flush after each write so results are recoverable if simulation crashes mid-run
    file_results.parameters["flush_output"] = True
    # all fields share the same mesh — avoids redundant mesh storage in the XDMF file
    file_results.parameters["functions_share_mesh"] = True

if(save_xdmf == 1):
    # write the undeformed initial state as the t=0 snapshot
    writeResults(BVP, file_results, t=0.0)

# =========================================================
# Storage array setup
# =========================================================
# totSteps includes step 0 (initial state) plus all loading steps
totSteps = steps + 1
# Forc matrix: rows = steps, columns = [disp, volume, strain_energy, inner_area, remote_p, inner_p]
Forc = np.zeros([totSteps, 6])
ii = 0

Forc[ii, 0] = disp_tot * (t / T_f)

# 2D deformation gradient and its cofactor for computing deformed cavity volume via surface integral
if(simulation_type != 1):
    F2D = as_tensor([[BVP.F_3D[0,0], BVP.F_3D[0,1]],
                     [BVP.F_3D[1,0], BVP.F_3D[1,1]]])
    J2D = det(F2D)
    # cofactor of F2D: used to map reference normals to deformed normals (Nanson's formula in 2D)
    F2Dcof = J2D * inv(F2D.T)
    n_dummy = dot(F2Dcof, BVP.Facet_normal)

# traction vector on the deformed surface: T · n_deformed
traction = dot(BVP.T, BVP.deformed_unit_normal)
# normal component of traction = pressure (positive = compression in solid mechanics convention)
Pressure = dot(traction, BVP.deformed_unit_normal)

# integrate normal force over the outer boundary (used to compute remote pressure)
OuterNormalForce = 2*pi*Pressure*BVP.ArealStretch*BVP.x[0]*BVP.dsm(Outer_boundary_surface)
# integrate normal force over the cavity surface (used to check traction-free BC)
InnerNormalForce = 2*pi*Pressure*BVP.ArealStretch*BVP.x[0]*BVP.dsm(Inner_boundary_surface)
# deformed cavity volume: computed from surface integral using the divergence theorem
Volume = pi*pow((BVP.x[0]+BVP.u[0]),2)*(-n_dummy[0])*BVP.dsm(Inner_boundary_surface)
# surface energy integrated over the deformed cavity surface
Energy_surface = 2*pi*BVP.psi_gamma*BVP.x[0]*BVP.dsm(Inner_boundary_surface)
# deformed area of the outer boundary
OuterArea = 2*pi*BVP.ArealStretch*BVP.x[0]*BVP.dsm(Outer_boundary_surface)
# deformed area of the cavity surface
InnerArea = 2*pi*BVP.ArealStretch*BVP.x[0]*BVP.dsm(Inner_boundary_surface)
# total bulk strain energy integrated over the domain
Energy_bulk = 2*pi*BVP.psi_bulk*BVP.x[0]*BVP.dxm

def V_cavity_theoretical(simulation_type, a_param):
    """
    Compute the theoretical cavity volume for a sphere given its radius a_param.

    Parameters
    ----------
    simulation_type : int
        2 = axisymmetric (full sphere volume = 4/3 * pi * a^3)
    a_param : float
        Current cavity radius
    """
    if simulation_type == 1:
        return (pi / 6) * a_param**3    # sphere quadrant (3D quarter model)
    elif simulation_type == 2:
        return (4 * pi / 3) * a_param**3  # full sphere (axisymmetric model)
    elif simulation_type == 3:
        return (pi / 4) * a_param**2    # quarter circle cross-section (plane strain)

# verify that the numerical initial cavity volume matches the analytical value
if(print_volume_check == 1):
    V0_theoretical = V_cavity_theoretical(simulation_type, A_param)
    V0_numerical = assemble(Volume)
    print(f"Initial cavity volume numerical: {V0_numerical:.6f}")
    print(f"Theoretical inital cavity volume: {V0_theoretical:.6f}")

# store initial state (t=0)
Forc[ii, 1] = assemble(Volume)
Forc[ii, 2] = assemble(Energy_bulk)
Forc[ii, 3] = assemble(InnerArea)
# remote pressure = total outer normal force / outer deformed area
Forc[ii, 4] = assemble(OuterNormalForce) / assemble(OuterArea)
# inner pressure = total cavity normal force / cavity deformed area (should stay ~0: traction-free)
Forc[ii, 5] = assemble(InnerNormalForce) / assemble(InnerArea)

# =========================================================
# Incremental loading loop
# =========================================================
while (t < T_f):

    ii += 1

    try:
        (iter, converged) = solver.solve()
        if not converged:
            # save last converged state before exiting if solver fails mid-run
            if save_xdmf == 1:
                writeResults(BVP, file_results, t - dt)
            raise RuntimeError("Solver did not converge at time {}".format(t))
    except RuntimeError as e:
        print(f"Solver failed at time {t}: {e}")
        break

    # write XDMF output every 1000 steps and at the final step
    if save_xdmf == 1 and (ii % 1000 == 0 or ii == steps):
        writeResults(BVP, file_results, t)

    # print progress every 10 steps
    if (ii-1) % 10 < 0.1:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        if(remote_displacement == 1):
            print("Step: Remote displacement |   Simulation Time: {}    |     Wallclock Time: {}".format(t, current_time))
        else:
            print("Step: Inner displacement   |   Simulation Time: {}    |     Wallclock Time: {}".format(t, current_time))

        print("Newton iterations: {}".format(iter))

        if(print_volume_check == 1):
            if(remote_displacement == 1):
                uy_theoretical = disp_tot * (t / T_f)
                # theoretical cavity radius from volume conservation: V = (4/3)*pi*a^3 = (4/3)*pi*A^3 + pi*B^2*uy
                a_theoretical = A_param * ((0.75*(BbyA**3.0)*uy_theoretical/B_param)+1)**(1./3.0)
            else:
                a_theoretical = A_param + disp_tot * (t / T_f)
                uy_theoretical = (4./3.) * B_param * (((a_theoretical**3.0/A_param**3.0)-1)/BbyA**3.0)
            V_theoretical = V_cavity_theoretical(simulation_type, a_theoretical)
            V_numerical = assemble(Volume)
            print(f"Cavity volume numerical: {V_numerical:.6f}")
            print(f"Theoretical cavity volume: {V_theoretical:.6f}")

            if(detailed_check == 1):
                print(f"Theoretical a: {a_theoretical}")
                print(f"Theoretical remote uy: {uy_theoretical}")
                # probe deformed coordinates at key points on the cavity boundary
                A_x = (A_param, y_c)             # equatorial point of cavity
                A_y_top = (0.0, y_c + A_param)   # top pole of cavity
                A_y_bottom = (0.0, y_c - A_param) # bottom pole of cavity (nearest to symmetry plane)
                B_1 = (0.0, B_param + y_c)        # top boundary on axis
                B_2 = (B_param, B_param + y_c)    # top boundary at right edge
                B_3 = (0.5*B_param, B_param + y_c) # top boundary at midpoint

                w_disp = BVP.w.sub(0)
                print(f"Deformed x at A_x: {A_param + w_disp.sub(0)(*A_x)}")
                print(f"Deformed y at A_y_top: {y_c + A_param + w_disp.sub(1)(*A_y_top)}")
                print(f"Deformed y at A_y_bottom: {y_c - A_param + w_disp.sub(1)(*A_y_bottom)}")
                print(f"ydisp at (0,B+y_c): {w_disp.sub(1)(*B_1)}")
                print(f"ysdisp at (B,B+y_c): {w_disp.sub(1)(*B_2)}")
                print(f"ysdisp at (0.5B,B+y_c): {w_disp.sub(1)(*B_3)}")

    # store results for this step
    Forc[ii, 0] = disp_tot * (t / T_f)
    Forc[ii, 1] = assemble(Volume)
    Forc[ii, 2] = assemble(Energy_bulk)
    Forc[ii, 3] = assemble(InnerArea)
    Forc[ii, 4] = assemble(OuterNormalForce) / assemble(OuterArea)
    Forc[ii, 5] = assemble(InnerNormalForce) / assemble(InnerArea)

    # update DOF vector for next step (store current solution as "old" state)
    BVP.w_old.vector()[:] = BVP.w.vector()

    # advance pseudo-time
    t += dt
    # update the time variable inside the time-dependent BC expressions
    if(remote_displacement == 0):
        dispx.t = t
    dispy.t = t

# =========================================================
# Final step summary
# =========================================================
if(print_volume_check == 1):
    if(remote_displacement == 1):
        uy_theoretical = disp_tot * (t / T_f)
        a_theoretical = A_param * ((0.75*(BbyA**3.0)*uy_theoretical/B_param)+1)**(1./3.0)
    else:
        a_theoretical = A_param + disp_tot * (t / T_f)
        uy_theoretical = (4./3.) * B_param * (((a_theoretical**3.0/A_param**3.0)-1)/BbyA**3.0)

    V_theoretical = V_cavity_theoretical(simulation_type, a_theoretical)
    V_numerical = assemble(Volume)

    print("\nFinal Step Results:")
    print("-" * 50)
    print("Step: Remote displacement |   Simulation Time: {}    |     Wallclock Time: {}".format(t, current_time))
    print("Newton iterations: {}".format(iter))
    print(f"Cavity volume numerical: {V_numerical:.6f}")
    print(f"Theoretical cavity volume: {V_theoretical:.6f}")

    if(detailed_check == 1):
        print(f"Theoretical a: {a_theoretical}")
        print(f"Theoretical remote uy: {uy_theoretical}")
        A_x = (A_param, y_c)
        A_y_top = (0.0, y_c + A_param)
        A_y_bottom = (0.0, y_c - A_param)
        B_1 = (0.0, B_param + y_c)
        B_2 = (B_param, B_param + y_c)
        B_3 = (0.5*B_param, B_param + y_c)

        w_disp = BVP.w.sub(0)
        print(f"Deformed x at A_x: {A_param + w_disp.sub(0)(*A_x)}")
        print(f"Deformed y at A_y_top: {y_c + A_param + w_disp.sub(1)(*A_y_top)}")
        print(f"Deformed y at A_y_bottom: {y_c - A_param + w_disp.sub(1)(*A_y_bottom)}")
        print(f"ydisp at (0,B+y_c): {w_disp.sub(1)(*B_1)}")
        print(f"ysdisp at (B,B+y_c): {w_disp.sub(1)(*B_2)}")
        print(f"ysdisp at (0.5B,B+y_c): {w_disp.sub(1)(*B_3)}")

print("End computation")

endTime = datetime.now()
elapseTime = endTime - startTime
print("------------------------------------")
print("Elapsed real time:  {}".format(elapseTime))
print("------------------------------------")

# =========================================================
# Save output data
# =========================================================
if(store_data_pickle == 1 or store_data_mat == 1):
    print("-----------------------------------------")
    print("Saving data")

    mdic = {"storage_var": Forc}
    if(store_data_mat == 1):
        print("Saving data in .mat file")
        savenamemat = savename + ".mat"
        # save Forc matrix as a .mat file: each row is a load step, columns as described in header
        scipy.io.savemat(savenamemat, mdic)

    if(store_data_pickle == 1):
        print("Saving data in .pickle file")
        savenamepickle = savename + ".pickle"
        with open(savenamepickle, 'wb') as f:
            pickle.dump(Forc, f)

    print("-----------------------------------------")