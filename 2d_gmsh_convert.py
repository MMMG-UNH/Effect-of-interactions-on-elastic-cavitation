#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert the mesh to .xdmf format using meshio

first step: gmsh -2 -format msh2 2Cavity_yc20.00_Finest.geo
second step: python3 2d_gmsh_convert.py 2Cavity_yc20.00_Finest

This script reads a GMSH .msh file and converts it to XDMF format for use in FEniCS.
It extracts the triangle and line (facet) cells from the mesh and saves them as separate XDMF files.
This script is intended for use with 2D meshes only.
"""

# pass filename without .msh extension as a command-line argument
import sys
import numpy
# the mesh filename (without extension) is passed as the first argument
filename = str(sys.argv[1])
# construct full file paths for input (.msh) and both output (.xdmf) files
filenamemsh = (filename + ".msh")
# main XDMF file: stores the triangular domain mesh used by FEniCS for the bulk solve
filenamexdmf = (filename + ".xdmf")
# facet XDMF file: stores the boundary line elements used by FEniCS for applying boundary conditions
filenamexdmf_facet = ("facet_" + filename + ".xdmf")

import meshio
# read the Gmsh .msh file into a meshio mesh object
mesh_from_file = meshio.read(filenamemsh)

# print mesh statistics to verify the mesh was read correctly before conversion
print("\nDetailed Mesh Information:")
print("-" * 50)
print(f"Points: {len(mesh_from_file.points)}")
print("\nCell Data:")
for key in mesh_from_file.cell_data_dict:
    print(f"\nKey: {key}")
    for cell_type, data in mesh_from_file.cell_data_dict[key].items():
        # unique values correspond to the physical group tags assigned in the .geo file
        print(f"  {cell_type}: unique values = {numpy.unique(data)}")


"""
Extract cells and boundary data.

Now that we have created the mesh, we need to extract the cells 
and physical data. We need to create a separate file for the 
facets (lines),  which we will use when we define boundary 
conditions in  Fenics. We do this  with the following convenience 
function. Note that as we would like a  2 dimensional mesh, we need to 
remove the z-values in the mesh coordinates, if any.
"""
import numpy
def create_mesh(mesh, cell_type, prune_z=False):
    # extract only cells of the specified type (e.g. "triangle" or "line")
    cells = mesh.get_cells_type(cell_type)
    # retain the physical group tag (boundary/region ID) for each cell
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    # drop the z-coordinate for 2D meshes to avoid FEniCS dimension mismatch
    points = mesh.points[:,:2] if prune_z else mesh.points
    # build a new minimal meshio Mesh containing only this cell type and its tags
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells},\
                           cell_data={"name_to_read":[cell_data]})
    return out_mesh

"""
With this function in hand, we can save the facet line mesh 
and the domain triangle  mesh in `XDMF` format 
"""

# extract and save boundary facets (line elements) as a separate XDMF file
# FEniCS reads this file to identify boundary regions for Dirichlet BCs
line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write(filenamexdmf_facet, line_mesh)

# extract and save the domain mesh (triangular elements) as the main XDMF file
# FEniCS reads this file to build the function spaces and assemble the weak form
triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write(filenamexdmf, triangle_mesh)