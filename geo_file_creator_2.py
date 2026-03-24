import numpy as np
from scipy.optimize import fsolve
import sys

# ==== USAGE & ARGS ====
# enforce that exactly two command-line arguments are provided (yc and mesh_quality)
if len(sys.argv) != 3:
    print("Usage: python3 geo_file_creator.py <yc> <mesh_quality>")
    print("  <mesh_quality>: Finest, Fine, Medium, or Coarse")
    sys.exit(1)

# yc: dimensionless cavity center distance from the bottom symmetry plane (d/R value)
yc = float(sys.argv[1])
# capitalize to match dictionary keys regardless of user input case
quality = sys.argv[2].capitalize()

# mapping mesh quality label to (h_min, res_med, res_low) — element sizes near cavity, mid-field, and far-field
quality_params = {
    'Finest': (0.002, 0.125, 0.5),
    'Fine'  : (0.005, 0.25,  1.0),
    'Medium': (0.01,  0.5,   2.0),
    'Coarse': (0.02,  1.0,   4.0)
}
# abort with a clear message if the provided quality label is not recognized
if quality not in quality_params:
    raise ValueError(f"Unknown mesh quality '{quality}'. Choose from {list(quality_params)}.")

# unpack the three resolution levels for use in .geo point definitions
resolution_high, resolution_med, resolution_low = quality_params[quality]
# h_min is the smallest element size, applied at and near the cavity surface
h_min = resolution_high

# fixed geometry parameters (in units of cavity radius R = a = 0.5)
# N_factor controls the number of transfinite mesh divisions per unit length along Line(3)
N_factor = 25
# L, H: half-width and height of the rectangular simulation domain (= 100R each)
L = 50.0
H = 50.0
# a: undeformed cavity radius
a = 0.5
# yc_crit: threshold above which a graded transfinite line is needed along the axis below the cavity
yc_crit = 2.5

# ==== BUILD .geo LINES ====
# initialize the list that will hold all Gmsh .geo script lines
geo_lines = []

# write domain dimensions and resolution parameters as Gmsh variables
geo_lines.extend([
    "// Domain dimensions",
    f"L = {L};",
    f"H = {H};",
    f"a = {a};  // cavity radius",
    f"yc = {yc}; // cavity center y-position\n",
    "// Mesh resolution (h_min around cavity, medium, coarse)",
    f"resolution_high = {resolution_high};",
    f"resolution_med  = {resolution_med};",
    f"resolution_low  = {resolution_low};\n"
])

# for large yc, the axis below the cavity (Line 3) spans a large distance;
# a transfinite graded progression is used so elements stay fine near the cavity
# and coarsen toward the bottom without excessive total element count
if yc >= yc_crit:
    # N: total number of transfinite divisions along Line(3)
    N = int(np.floor(N_factor * yc))
    # solve for the geometric progression ratio r such that the first element = h_min
    # and the total length of N elements sums to yc
    def f(r): return h_min * (1 - r**N) / (1 - r) - yc
    r_solution = fsolve(f, 1.1)[0]
    geo_lines += [
        "// Derived mesh control",
        f"N = {N};",
        f"r = {r_solution:.6f};\n"
    ]

# define the 7 key geometric points of the axisymmetric domain
# all points lie on the symmetry axis (x=0) or domain boundary
geo_lines += [
    "// Points",
    "Point(1) = {0, H + yc, 0, resolution_low};         // top-left",
    "Point(2) = {0, yc + a, 0, resolution_high};         // cavity top",
    "Point(3) = {0, yc,       0, resolution_high};       // cavity center",
    "Point(4) = {0, yc - a,   0, resolution_high};       // cavity bottom",
    f"Point(5) = {{0, 0, 0, {resolution_med}}};          // bottom-left",
    "Point(6) = {L, 0, 0, resolution_low};               // bottom-right",
    "Point(7) = {L, H + yc, 0, resolution_low};         // top-right\n"
]

# define boundary lines and the cavity arc
geo_lines += [
    "// Lines",
    "Line(1) = {1, 2};              // left wall above cavity",
    # Circle arc defined by: start point, center point, end point
    "Circle(2) = {4, 3, 2};         // cavity arc",
    "Line(3) = {4, 5};              // from cavity bottom to bottom-left"
]
# apply transfinite meshing with geometric progression to Line(3) only when yc is large
if yc >= yc_crit:
    geo_lines.append("Transfinite Line {3} = N Using Progression r;")
geo_lines += [
    "Line(4) = {5, 6};              // bottom edge",
    "Line(5) = {6, 7};              // right edge",
    "Line(6) = {7, 1};              // top edge\n"
]

# define the closed curve loop (note: -2 reverses the cavity arc orientation for consistent winding)
# then define the planar surface enclosed by this loop
geo_lines += [
    "// Curve Loop and Surface",
    "Curve Loop(8) = {1, -2, 3, 4, 5, 6};",
    "Plane Surface(9) = {8};\n",
    # physical groups assign boundary tags used later in FEniCS for applying boundary conditions
    "// Physical groups",
    'Physical Curve("left-top", 10)   = {1};',
    'Physical Curve("cavity", 11)     = {2};',
    'Physical Curve("left-bottom", 12) = {3};',
    'Physical Curve("bottom", 13)     = {4};',
    'Physical Curve("right", 14)      = {5};',
    'Physical Curve("top", 15)        = {6};',
    'Physical Surface("Body", 16)     = {9};'
]

# write all collected .geo lines to a file named after the yc value and mesh quality
filename = f"2Cavity_yc{yc:.2f}_{quality}.geo"
with open(filename, "w") as f:
    f.write("\n".join(geo_lines))

print(f"Generated '{filename}' with mesh quality '{quality}' and yc={yc:.2f}.")